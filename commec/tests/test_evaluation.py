#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Data-driven, end-to-end detection ("evaluation") tests for ``commec screen``.

Unlike the mocked logic tests (see ``screen_factory.py``), these tests run the *real*
screening pipeline against *real* reference databases and assert on the detection outcome
for known sequences. Each test case is a subfolder of ``test_data/eval/`` containing:

  * a ``metadata.yaml`` (or ``.yml``) file describing the case and its expected outcomes;
  * one FASTA file of input sequence(s) to screen.

A FASTA is expected to hold a single record by default, but may hold several (for example
to test "distribution attacks" that split a sequence of concern across multiple records);
expected outcomes are keyed by FASTA record id, so multi-record cases are supported.

See ``test_data/eval/README.md`` for the full ``metadata.yaml`` schema.

Because the real reference databases are not stored in the repository, these tests are
SKIPPED unless a database directory is supplied, either via::

    pytest --commec-db-dir /path/to/commec-dbs

or by setting the ``COMMEC_DB_DIR`` environment variable.
"""

import os
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from commec.config.json_io import get_screen_data_from_json
from commec.config.result import QueryResult, ScreenResult, ScreenStatus
from commec.screen import ScreenArgumentParser, add_args, run

# Root directory under which each subfolder is one evaluation case.
EVAL_DIR = Path(__file__).parent / "test_data" / "eval"

# Accepted names for a case's metadata file, in order of preference.
METADATA_FILENAMES = ("metadata.yaml", "metadata.yml")

# Per-step status fields on QueryScreenStatus that a case may assert against.
# "screen_status" is the overall outcome and is the only required expectation.
STATUS_FIELDS = (
    "screen_status",
    "biorisk",
    "protein_taxonomy",
    "nucleotide_taxonomy",
    "low_concern",
)


# --------------------------------------------------------------------------------------
# Case discovery
# --------------------------------------------------------------------------------------
def _discover_cases() -> list[Path]:
    """
    Return the sorted list of evaluation case directories under ``EVAL_DIR``.

    A directory is a case if it contains a metadata file. Folders whose name starts with
    ``_`` or ``.`` are ignored, so e.g. ``_template/`` can hold documentation/examples
    without being collected as a test.
    """
    if not EVAL_DIR.is_dir():
        return []

    cases = []
    for entry in sorted(EVAL_DIR.iterdir()):
        if not entry.is_dir() or entry.name.startswith(("_", ".")):
            continue
        if any((entry / name).is_file() for name in METADATA_FILENAMES):
            cases.append(entry)
    return cases


def _case_id(case_dir: Path) -> str:
    return case_dir.name


_CASES = _discover_cases()

if _CASES:
    _CASE_PARAMS = [pytest.param(c, id=_case_id(c)) for c in _CASES]
else:
    # Keep a single, clearly-labelled skipped test rather than collecting nothing, so
    # that an empty/absent eval directory is visible in the test report.
    _CASE_PARAMS = [
        pytest.param(
            None,
            id="no-eval-cases",
            marks=pytest.mark.skip(reason=f"No evaluation cases found in {EVAL_DIR}"),
        )
    ]


# --------------------------------------------------------------------------------------
# Fixtures and helpers
# --------------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def eval_database_dir(request) -> Path:
    """
    Resolve the real commec database directory from ``--commec-db-dir`` or the
    ``COMMEC_DB_DIR`` env var. Skips the test when neither is set or the path is invalid.
    """
    raw = request.config.getoption("--commec-db-dir") or os.environ.get("COMMEC_DB_DIR")
    if not raw:
        pytest.skip(
            "No commec database directory provided; pass --commec-db-dir or set"
            " COMMEC_DB_DIR to run the end-to-end detection tests."
        )
    db_dir = Path(raw).expanduser()
    if not db_dir.is_dir():
        pytest.skip(f"commec database directory does not exist: {db_dir}")
    return db_dir


def _load_metadata(case_dir: Path) -> dict:
    """Load and minimally validate the metadata for a case directory."""
    meta_path = next(
        (case_dir / name for name in METADATA_FILENAMES if (case_dir / name).is_file()),
        None,
    )
    assert meta_path is not None, f"No metadata file found in {case_dir}"

    with open(meta_path, encoding="utf-8") as handle:
        metadata = yaml.safe_load(handle) or {}

    assert isinstance(metadata, dict), f"{meta_path} must contain a YAML mapping"
    assert "expected" in metadata, f"{meta_path} is missing the required 'expected' key"
    return metadata


def _resolve_fasta(case_dir: Path, metadata: dict) -> Path:
    """
    Determine which FASTA in the case directory to screen.

    Uses ``metadata['fasta']`` if given; otherwise requires exactly one ``*.fasta`` /
    ``*.fa`` file in the directory.
    """
    named = metadata.get("fasta")
    if named:
        fasta = case_dir / named
        assert fasta.is_file(), f"FASTA named in metadata not found: {fasta}"
        return fasta

    candidates = sorted(
        p for p in case_dir.iterdir() if p.suffix.lower() in (".fasta", ".fa")
    )
    assert candidates, f"No FASTA file found in {case_dir}"
    assert len(candidates) == 1, (
        f"Multiple FASTA files in {case_dir}; name one explicitly with the 'fasta' key"
        f" in metadata. Found: {[p.name for p in candidates]}"
    )
    return candidates[0]


def _parse_status(value) -> ScreenStatus:
    """
    Coerce a metadata status value into a ``ScreenStatus``.

    Accepts either the enum *value* (the human string, e.g. ``"Flag"``,
    ``"Flag (Cleared)"``) or the enum *name* (e.g. ``"FLAG"``, ``"CLEARED_FLAG"``).
    """
    if isinstance(value, ScreenStatus):
        return value
    text = str(value)
    try:
        return ScreenStatus(text)
    except ValueError:
        pass
    try:
        return ScreenStatus[text.upper()]
    except KeyError as exc:
        valid = ", ".join(f"{s.name}={s.value!r}" for s in ScreenStatus)
        raise AssertionError(
            f"Unknown ScreenStatus {value!r} in metadata. Valid options: {valid}"
        ) from exc


def _normalise_expectations(expected) -> dict[str, dict]:
    """
    Normalise the ``expected`` block into ``{query_id: {field: expectation, ...}}``.

    Supports a single-query shorthand where the status fields are written directly under
    ``expected`` (with no query-id nesting); these are applied to the sole query under the
    sentinel key ``None``.
    """
    assert isinstance(expected, dict), "'expected' must be a mapping"

    # Shorthand: status fields written directly (single, unnamed query).
    if set(expected).issubset(set(STATUS_FIELDS) | {"hits", "rationale_contains"}):
        return {None: expected}

    # Otherwise, keys are query ids mapping to their own expectation blocks.
    out = {}
    for query_id, fields in expected.items():
        assert isinstance(fields, dict), (
            f"Expectation for query {query_id!r} must be a mapping of fields"
        )
        out[str(query_id)] = fields
    return out


def _find_query(result: ScreenResult, query_id) -> QueryResult:
    """
    Locate the QueryResult for an expected query id, tolerant of how commec keys queries.

    ``query_id`` may be ``None`` (the single-query shorthand), the internal query name
    (the ``queries`` dict key), or the original FASTA record id.
    """
    queries = result.queries
    if query_id is None:
        assert len(queries) == 1, (
            "Single-query shorthand used, but the screen produced"
            f" {len(queries)} queries: {list(queries)}. Key 'expected' by query id."
        )
        return next(iter(queries.values()))

    if query_id in queries:
        return queries[query_id]

    for query in queries.values():
        if query.query == query_id:
            return query

    resolved = result.get_query(query_id)
    assert resolved is not None, (
        f"Expected query {query_id!r} not found in screen output."
        f" Queries present: {list(queries)}"
    )
    return resolved


def _run_screen(
    fasta: Path, db_dir: Path, output_dir: Path, screen_args: list[str]
) -> ScreenResult:
    """Run a real ``commec screen`` and return the parsed ScreenResult."""
    arguments = [
        "commec-screen",
        str(fasta),
        "-d",
        str(db_dir),
        "-o",
        str(output_dir),
        "-F",  # force a fresh run into the (empty) tmp output dir
    ]
    arguments.extend(screen_args)

    with patch("sys.argv", arguments):
        parser = ScreenArgumentParser()
        add_args(parser)
        args = parser.parse_args()
        run(args)

    # When -o is a directory, the output JSON is named after the input FASTA.
    output_json = output_dir / f"{fasta.stem}.output.json"
    if not output_json.is_file():
        # Fall back to any *.output.json produced, to be robust to naming changes.
        produced = sorted(output_dir.glob("*.output.json"))
        assert produced, f"No screen JSON output produced in {output_dir}"
        output_json = produced[0]

    return get_screen_data_from_json(str(output_json))


# --------------------------------------------------------------------------------------
# The test
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("case_dir", _CASE_PARAMS)
def test_evaluation_case(case_dir: Path, eval_database_dir: Path, tmp_path: Path):
    """
    Run ``commec screen`` on a single evaluation case and assert its expected outcomes.
    """
    metadata = _load_metadata(case_dir)
    fasta = _resolve_fasta(case_dir, metadata)
    screen_args = [str(a) for a in metadata.get("screen_args", [])]

    result = _run_screen(fasta, eval_database_dir, tmp_path, screen_args)

    expectations = _normalise_expectations(metadata["expected"])
    case_name = case_dir.name

    failures: list[str] = []
    for query_id, fields in expectations.items():
        query = _find_query(result, query_id)
        label = query.query or query_id or "<single query>"

        for field in STATUS_FIELDS:
            if field not in fields:
                continue
            expected_status = _parse_status(fields[field])
            actual_status = getattr(query.status, field)
            if actual_status != expected_status:
                failures.append(
                    f"[{case_name}:{label}] {field}: expected {expected_status.value!r},"
                    f" got {actual_status.value!r}"
                )

        # Optional: assert that named hits were detected for this query.
        for hit_name in fields.get("hits", []):
            if hit_name not in query.hits:
                failures.append(
                    f"[{case_name}:{label}] expected hit {hit_name!r} not found."
                    f" Hits present: {list(query.hits)}"
                )

        # Optional: assert substrings appear in the rationale text.
        for fragment in fields.get("rationale_contains", []):
            if fragment not in query.status.rationale:
                failures.append(
                    f"[{case_name}:{label}] rationale missing {fragment!r}."
                    f" Rationale: {query.status.rationale!r}"
                )

    assert not failures, "Evaluation case mismatch:\n" + "\n".join(failures)
