#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Tests for commec/screeners/check_reg_path.py

Covers each refactored function in logical order: input validation, data
preparation, region annotation collection, enrichment, outcome determination,
hit orchestration, ScreenResult mutation, and the main entry point.
"""
import pytest
from unittest.mock import MagicMock, patch, call
import pandas as pd

from commec.screeners.check_reg_path import (
    _check_inputs,
    _mark_initial_hits,
    _collect_region_annotations,
    _extract_unique_annotations,
    _enrich_regulated_annotations,
    _enrich_non_regulated_annotations,
    _determine_screen_outcome,
    _process_regulated_hit,
    _apply_hit_to_screen_result,
    _update_using_control_lists,
    parse_taxonomy_hits,
)
from commec.config.result import (
    ScreenResult,
    ScreenStatus,
    ScreenStep,
    HitResult,
    HitScreenStatus,
    MatchRange,
    TaxonomyAnnotation,
    QueryResult,
    QueryScreenStatus,
)
from commec.control_list.containers import ControlListOutput, ListMode


# ── Helpers ────────────────────────────────────────────────────────────────────

def _blast_row(**overrides) -> dict:
    """Return a minimal BLAST-output row dict suitable for building test DataFrames."""
    row = {
        "query acc.": "seq1",
        "subject title": "Dangerous virus polymerase",
        "subject acc.": "NC_001234",
        "subject tax ids": 11234,
        "evalue": 1e-10,
        "bit score": 500.0,
        "% identity": 99.0,
        "query length": 1000,
        "q. start": 100,
        "q. end": 300,
        "subject length": 300,
        "s. start": 1,
        "s. end": 200,
        "regulated": True,
        "log evalue": 10.0,
        "q. coverage": 0.2,
        "s. coverage": 1.0,
        "species": "",
        "genus": "",
        "category": "",
        "list_acronym": "",
    }
    row.update(overrides)
    return row


def _make_df(*rows) -> pd.DataFrame:
    """Build a DataFrame from one or more row dicts from _blast_row()."""
    return pd.DataFrame([_blast_row(**r) for r in rows])


def _mock_control_list(status: ListMode = ListMode.COMPLIANCE, acronym: str = "DTL") -> MagicMock:
    cl = MagicMock()
    cl.status = status
    cl.acronym = acronym
    return cl


def _mock_control_output(category: str = "Viruses", species: str = "Virus sp.",
                          genus: str = "Virusgenus", list_str: str = "DTL") -> ControlListOutput:
    out = ControlListOutput()
    out.category = category
    out.species = species
    out.genus = genus
    out.list = list_str
    return out


def _make_screen_result(query_name: str = "seq1") -> ScreenResult:
    result = ScreenResult()
    result.queries[query_name] = QueryResult(query=query_name)
    return result


# ── _check_inputs ──────────────────────────────────────────────────────────────

class TestCheckInputs:

    def _handler(self, validates: bool) -> MagicMock:
        h = MagicMock()
        h.validate_output.return_value = validates
        h.out_file = "/mock/out.blast"
        return h

    def test_returns_false_when_handler_output_invalid(self):
        assert _check_inputs(self._handler(False), "/lc.csv", "/br.csv", "/taxdir") is False

    def test_returns_false_when_low_concern_missing(self, tmp_path):
        taxdir = tmp_path / "taxdir"
        taxdir.mkdir()
        assert _check_inputs(self._handler(True), "/nonexistent.csv", "/br.csv", str(taxdir)) is False

    def test_returns_false_when_taxonomy_dir_missing(self, tmp_path):
        lc = tmp_path / "lc.csv"
        lc.write_text("1234\n")
        assert _check_inputs(self._handler(True), str(lc), "/br.csv", "/nonexistent_dir") is False

    def test_returns_true_when_all_valid(self, tmp_path):
        lc = tmp_path / "lc.csv"
        lc.write_text("1234\n")
        taxdir = tmp_path / "taxdir"
        taxdir.mkdir()
        assert _check_inputs(self._handler(True), str(lc), "/br.csv", str(taxdir)) is True


# ── _collect_region_annotations ───────────────────────────────────────────────

class TestCollectRegionAnnotations:

    def _make_region_and_query_data(self, q_start=100, q_end=300):
        """Create a simple region row and a surrounding query DataFrame."""
        region = pd.Series(_blast_row(**{"q. start": q_start, "q. end": q_end, "regulated": True}))
        query_data = pd.DataFrame([
            _blast_row(**{"q. start": q_start, "q. end": q_end, "regulated": True,
                          "subject acc.": "NC_001", "subject tax ids": 111}),
            _blast_row(**{"q. start": q_start, "q. end": q_end, "regulated": False,
                          "subject acc.": "NC_002", "subject tax ids": 222}),
            _blast_row(**{"q. start": 999, "q. end": 1999, "regulated": False,
                          "subject acc.": "NC_003", "subject tax ids": 333}),
        ])
        return region, query_data

    @patch("commec.screeners.check_reg_path._update_using_control_lists",
           side_effect=lambda r, nr: (r, nr))
    def test_shared_start_matched(self, _mock_ucl):
        region, query_data = self._make_region_and_query_data(q_start=100, q_end=300)
        reg, non_reg = _collect_region_annotations(region, query_data)
        # NC_001 is regulated at same coords, NC_003 has different coords
        assert 111 in reg["subject tax ids"].values
        assert 333 not in non_reg["subject tax ids"].values

    @patch("commec.screeners.check_reg_path._update_using_control_lists",
           side_effect=lambda r, nr: (r, nr))
    def test_non_regulated_capped_at_10(self, _mock_ucl):
        rows = [
            _blast_row(**{"q. start": 100, "q. end": 300, "regulated": False,
                          "subject tax ids": i, "subject acc.": f"NC_{i:04d}"})
            for i in range(20)
        ]
        region = pd.Series(_blast_row(**{"q. start": 100, "q. end": 300}))
        query_data = pd.DataFrame(rows)
        _reg, non_reg = _collect_region_annotations(region, query_data)
        assert len(non_reg) <= 10

    @patch("commec.screeners.check_reg_path._update_using_control_lists",
           side_effect=lambda r, nr: (r, nr))
    def test_deduplication_by_taxid(self, _mock_ucl):
        rows = [
            _blast_row(**{"q. start": 100, "q. end": 300, "regulated": True,
                          "subject tax ids": 111, "subject acc.": "NC_001"}),
            _blast_row(**{"q. start": 100, "q. end": 300, "regulated": True,
                          "subject tax ids": 111, "subject acc.": "NC_001b"}),  # same taxid
        ]
        region = pd.Series(_blast_row(**{"q. start": 100, "q. end": 300}))
        query_data = pd.DataFrame(rows)
        reg, _non_reg = _collect_region_annotations(region, query_data)
        # Only one row should survive deduplication on taxid
        assert len(reg) == 1


# ── _extract_unique_annotations ───────────────────────────────────────────────

class TestExtractUniqueAnnotations:

    def test_returns_correct_annotation_objects(self):
        df = _make_df(
            {"evalue": 1e-5, "subject tax ids": 999, "subject acc.": "NC_999",
             "subject title": "Test hit"},
        )
        result = _extract_unique_annotations(df)
        assert len(result) == 1
        ann = next(iter(result))
        assert isinstance(ann, TaxonomyAnnotation)
        assert ann.taxid == 999
        assert ann.target_hit == "NC_999"

    def test_deduplicates_identical_rows(self):
        df = pd.DataFrame([
            _blast_row(**{"evalue": 1e-5, "subject tax ids": 111}),
            _blast_row(**{"evalue": 1e-5, "subject tax ids": 111}),
        ])
        assert len(_extract_unique_annotations(df)) == 1

    def test_two_distinct_annotations(self):
        df = pd.DataFrame([
            _blast_row(**{"evalue": 1e-5, "subject tax ids": 111, "subject acc.": "NC_A"}),
            _blast_row(**{"evalue": 1e-5, "subject tax ids": 222, "subject acc.": "NC_B"}),
        ])
        assert len(_extract_unique_annotations(df)) == 2


# ── _enrich_regulated_annotations ─────────────────────────────────────────────

class TestEnrichRegulatedAnnotations:

    def _annotation(self, taxid: int = 111) -> dict:
        return {"taxid": str(taxid), "evalue": 1e-5,
                "target_hit": "NC_001", "target_description": "Test"}

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_populates_control_list_field(self, mock_gcl, mock_gr):
        cl_out = _mock_control_output(category="Viruses")
        mock_gr.return_value = ([cl_out], [])
        mock_gcl.return_value = _mock_control_list()
        annotations = [self._annotation()]
        result, _ = _enrich_regulated_annotations(annotations)
        assert "control_list" in result[0]

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_collects_domain_and_breaks_after_first_major_domain(self, mock_gcl, mock_gr):
        # Two control entries: "Proteins" then "Viruses" — should break after "Viruses"
        cl_out_proteins = _mock_control_output(category="Proteins")
        cl_out_viruses = _mock_control_output(category="Viruses")
        cl_out_bacteria = _mock_control_output(category="Bacteria")
        mock_gr.return_value = ([cl_out_proteins, cl_out_viruses, cl_out_bacteria], [])
        mock_gcl.return_value = _mock_control_list()
        annotations = [self._annotation()]
        _, domains = _enrich_regulated_annotations(annotations)
        # "Proteins" and "Viruses" should be in domains, but NOT "Bacteria"
        assert "Proteins" in domains
        assert "Viruses" in domains
        assert "Bacteria" not in domains

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_domain_count_per_annotation(self, mock_gcl, mock_gr):
        cl_out = _mock_control_output(category="Bacteria")
        mock_gr.return_value = ([cl_out], [])
        mock_gcl.return_value = _mock_control_list()
        annotations = [self._annotation(111), self._annotation(222)]
        _, domains = _enrich_regulated_annotations(annotations)
        assert domains.count("Bacteria") == 2


# ── _determine_screen_outcome ──────────────────────────────────────────────────

class TestDetermineScreenOutcome:

    def test_flag_when_only_regulated(self):
        status, desc, _ = _determine_screen_outcome(
            [{"taxid": "111"}], [], ["Viruses"], "Ebola polymerase"
        )
        assert status == ScreenStatus.FLAG
        assert "Ebola polymerase" in desc

    def test_pass_when_mixed(self):
        status, desc, _ = _determine_screen_outcome(
            [{"taxid": "111"}], [{"taxid": "222"}], ["Viruses"], "Ebola polymerase"
        )
        assert status == ScreenStatus.PASS
        assert "Mix of" in desc

    def test_pass_when_no_regulated_annotations(self):
        status, desc, _ = _determine_screen_outcome([], [{"taxid": "222"}], [], "Any hit")
        assert status == ScreenStatus.PASS
        assert "Regionally controlled" in desc

    def test_domains_text_defaults_to_sequence_when_empty(self):
        _, _, domains_text = _determine_screen_outcome([], [], [], "hit")
        assert domains_text == "sequence"

    def test_domains_text_uses_unique_domains(self):
        _, _, domains_text = _determine_screen_outcome(
            [{"taxid": "111"}], [], ["Viruses", "Viruses", "Bacteria"], "hit"
        )
        # Should be deduplicated
        parts = set(domains_text.split(", "))
        assert parts == {"Viruses", "Bacteria"}


# ── _update_using_control_lists ─────────────────────────────────────────────────

class TestUpdateUsingControlLists:

    def _make_regulated(self, taxid: int = 111, count: int = 1) -> pd.DataFrame:
        rows = [_blast_row(**{"subject tax ids": taxid + i, "regulated": True}) for i in range(count)]
        return pd.DataFrame(rows)

    def _make_non_regulated(self) -> pd.DataFrame:
        return pd.DataFrame([_blast_row(**{"subject tax ids": 9999, "regulated": False})])

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_compliance_entry_stays_regulated(self, mock_gcl, mock_gr):
        cl_out = _mock_control_output()
        mock_gr.return_value = ([cl_out], [])
        mock_gcl.return_value = _mock_control_list(ListMode.COMPLIANCE)
        reg = self._make_regulated()
        non_reg = self._make_non_regulated()
        new_reg, new_non_reg = _update_using_control_lists(reg, non_reg)
        assert len(new_reg) == 1
        assert len(new_non_reg) == 1  # unchanged

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_single_conditional_moves_to_non_regulated(self, mock_gcl, mock_gr):
        cl_out = _mock_control_output()
        mock_gr.return_value = ([cl_out], [])
        mock_gcl.return_value = _mock_control_list(ListMode.CONDITIONAL_NUM)
        reg = self._make_regulated(count=1)
        non_reg = pd.DataFrame()
        new_reg, new_non_reg = _update_using_control_lists(reg, non_reg)
        # Single conditional: below threshold, gets moved
        assert len(new_reg) == 0
        assert len(new_non_reg) == 1

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_two_conditionals_stay_regulated(self, mock_gcl, mock_gr):
        """When CONDITIONAL_NUMBER_TO_ALLOW (2) is reached, treat all as regulated."""
        cl_out = _mock_control_output()
        mock_gr.return_value = ([cl_out], [])
        mock_gcl.return_value = _mock_control_list(ListMode.CONDITIONAL_NUM)
        reg = self._make_regulated(count=2)
        non_reg = pd.DataFrame()
        new_reg, new_non_reg = _update_using_control_lists(reg, non_reg)
        # Threshold met: nothing moves
        assert len(new_reg) == 2
        assert len(new_non_reg) == 0

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_species_and_category_populated(self, mock_gcl, mock_gr):
        cl_out = _mock_control_output(category="Bacteria", species="E. coli")
        mock_gr.return_value = ([cl_out], [])
        mock_gcl.return_value = _mock_control_list(ListMode.COMPLIANCE)
        reg = self._make_regulated()
        new_reg, _ = _update_using_control_lists(reg, pd.DataFrame())
        assert new_reg.iloc[0]["species"] == "E. coli"
        assert new_reg.iloc[0]["category"] == "Bacteria"


# ── _apply_hit_to_screen_result ───────────────────────────────────────────────

class TestApplyHitToScreenResult:

    def _make_query_result(self) -> QueryResult:
        return QueryResult(query="seq1", status=QueryScreenStatus())

    def _make_hit(self, name: str = "Dangerous virus polymerase",
                  status: ScreenStatus = ScreenStatus.FLAG) -> HitResult:
        return HitResult(
            recommendation=HitScreenStatus(status, ScreenStep.TAXONOMY_AA),
            name=name,
            description="test desc",
            ranges=[MatchRange(1e-10, 1, 200, 100, 300)],
            annotations={"domain": ["Viruses"], "regulated_taxonomy": [{}]},
        )

    def test_step_status_updated(self):
        qr = self._make_query_result()
        hit = self._make_hit(status=ScreenStatus.FLAG)
        _apply_hit_to_screen_result(
            qr, "NC_001234", hit, hit.ranges, ["Viruses"], {}, ScreenStatus.FLAG, ScreenStep.TAXONOMY_AA
        )
        assert qr.status.protein_taxonomy == ScreenStatus.FLAG

    def test_new_hit_added_to_query_result(self):
        qr = self._make_query_result()
        hit = self._make_hit()
        _apply_hit_to_screen_result(
            qr, "NC_001234", hit, hit.ranges, ["Viruses"], {}, ScreenStatus.FLAG, ScreenStep.TAXONOMY_AA
        )
        assert hit.name in qr.hits

    def test_duplicate_hit_does_not_create_second_entry(self):
        qr = self._make_query_result()
        hit = self._make_hit()
        _apply_hit_to_screen_result(
            qr, "NC_001234", hit, hit.ranges, ["Viruses"], {}, ScreenStatus.FLAG, ScreenStep.TAXONOMY_AA
        )
        _apply_hit_to_screen_result(
            qr, "NC_001234", hit, hit.ranges, ["Viruses"], {}, ScreenStatus.FLAG, ScreenStep.TAXONOMY_AA
        )
        assert len(qr.hits) == 1


# ── parse_taxonomy_hits (integration) ─────────────────────────────────────────

class TestParseTaxonomyHits:
    """
    Integration tests for the main entry point.
    Filesystem and external-DB calls are patched; the BLAST pipeline is
    exercised through a fabricated top_hits DataFrame.
    """

    STEP = ScreenStep.TAXONOMY_AA

    def _make_handler(self, has_hits: bool = True) -> MagicMock:
        h = MagicMock()
        h.validate_output.return_value = True
        h.has_hits.return_value = has_hits
        h.out_file = "/mock/out.blast"
        return h

    def _base_patches(self):
        """Context managers that eliminate all filesystem and DB dependencies."""
        return [
            patch("os.path.exists", return_value=True),
            patch("commec.screeners.check_reg_path.read_blast"),
            patch("commec.screeners.check_reg_path.get_controlled_labels"),
            patch("commec.screeners.check_reg_path.get_top_hits"),
        ]

    def test_returns_1_on_invalid_inputs(self, tmp_path):
        handler = MagicMock()
        handler.validate_output.return_value = False
        handler.out_file = "/mock/out"
        result = ScreenResult()
        with patch("os.path.exists", return_value=True):
            ret = parse_taxonomy_hits(
                handler, "/lc.csv", "/br.csv", "/taxdir",
                result, {}, self.STEP, 1,
            )
        assert ret == 1

    def test_returns_0_and_pass_when_no_hits(self, tmp_path):
        handler = self._make_handler(has_hits=False)
        data = _make_screen_result("seq1")
        with patch("os.path.exists", return_value=True):
            ret = parse_taxonomy_hits(
                handler, "/lc.csv", "/br.csv", "/taxdir",
                data, {}, self.STEP, 1,
            )
        assert ret == 0
        assert data.queries["seq1"].status.protein_taxonomy == ScreenStatus.PASS

    def test_returns_0_and_pass_when_no_regulated_hits(self):
        handler = self._make_handler(has_hits=True)
        data = _make_screen_result("seq1")
        blast_df = _make_df({"regulated": False})
        top_hits_df = _make_df({"regulated": False})

        with (
            patch("os.path.exists", return_value=True),
            patch("commec.screeners.check_reg_path.read_blast", return_value=blast_df),
            patch("commec.screeners.check_reg_path.get_controlled_labels", return_value=blast_df),
            patch("commec.screeners.check_reg_path.get_top_hits", return_value=top_hits_df),
        ):
            ret = parse_taxonomy_hits(
                handler, "/lc.csv", "/br.csv", "/taxdir",
                data, {}, self.STEP, 1,
            )
        assert ret == 0
        assert data.queries["seq1"].status.protein_taxonomy == ScreenStatus.PASS

    def test_flags_query_with_regulated_hit(self):
        handler = self._make_handler(has_hits=True)
        data = _make_screen_result("seq1")

        blast_df = _make_df({"query acc.": "seq1", "regulated": True, "subject tax ids": 11234})
        top_hits_df = blast_df.copy()

        mock_query = MagicMock()
        mock_query.nc_to_nt_query_coords = lambda x: x

        control_out = _mock_control_output(category="Viruses", species="Test virus sp.")
        mock_cl = _mock_control_list(ListMode.COMPLIANCE)

        with (
            patch("os.path.exists", return_value=True),
            patch("commec.screeners.check_reg_path.read_blast", return_value=blast_df),
            patch("commec.screeners.check_reg_path.get_controlled_labels", return_value=blast_df),
            patch("commec.screeners.check_reg_path.get_top_hits", return_value=top_hits_df),
            patch("commec.screeners.check_reg_path.get_regulation",
                  return_value=([control_out], [])),
            patch("commec.screeners.check_reg_path.get_control_lists",
                  return_value=mock_cl),
        ):
            ret = parse_taxonomy_hits(
                handler, "/lc.csv", "/br.csv", "/taxdir",
                data, {"seq1": mock_query}, ScreenStep.TAXONOMY_AA, 1,
            )

        assert ret == 0
        assert data.queries["seq1"].status.protein_taxonomy == ScreenStatus.FLAG

    def test_pass_for_mixed_regulated_and_non_regulated_hit(self):
        handler = self._make_handler(has_hits=True)
        data = _make_screen_result("seq1")

        # Two hits at same position: one regulated, one not
        blast_df = pd.DataFrame([
            _blast_row(**{"query acc.": "seq1", "regulated": True,
                          "subject tax ids": 111, "subject acc.": "NC_REG",
                          "subject title": "Regulated virus"}),
            _blast_row(**{"query acc.": "seq1", "regulated": False,
                          "subject tax ids": 222, "subject acc.": "NC_NOT",
                          "subject title": "Harmless bacterium"}),
        ])
        top_hits_df = blast_df.copy()

        mock_query = MagicMock()
        mock_query.nc_to_nt_query_coords = lambda x: x

        control_out = _mock_control_output(category="Viruses")
        mock_cl = _mock_control_list(ListMode.COMPLIANCE)

        with (
            patch("os.path.exists", return_value=True),
            patch("commec.screeners.check_reg_path.read_blast", return_value=blast_df),
            patch("commec.screeners.check_reg_path.get_controlled_labels", return_value=blast_df),
            patch("commec.screeners.check_reg_path.get_top_hits", return_value=top_hits_df),
            patch("commec.screeners.check_reg_path.get_regulation",
                  return_value=([control_out], [])),
            patch("commec.screeners.check_reg_path.get_control_lists",
                  return_value=mock_cl),
        ):
            ret = parse_taxonomy_hits(
                handler, "/lc.csv", "/br.csv", "/taxdir",
                data, {"seq1": mock_query}, ScreenStep.TAXONOMY_AA, 1,
            )

        assert ret == 0
        assert data.queries["seq1"].status.protein_taxonomy == ScreenStatus.PASS
