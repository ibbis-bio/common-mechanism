#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Tests for commec/screeners/check_reg_path.py

Covers the functions that remain after the module refactor, in logical order:
input validation, hit-info formatting, annotation-level outcome determination,
cluster-level hit creation, per-query aggregation, and the main entry point.
"""
import pytest
from dataclasses import asdict
import json
import logging
from unittest.mock import MagicMock, patch
import pandas as pd

from commec.screeners.check_reg_path import (
    _create_hit_info,
    _create_hit_result_from_annotations,
    _create_hit_result_for_cluster,
    _get_hit_result_from_data,
    _build_log_message,
    parse_taxonomy_hits,
)
from commec.config.result import (
    ScreenResult,
    ScreenStatus,
    ScreenStep,
    HitResult,
    HitScreenStatus,
    MatchRange,
    QueryResult,
    QueryScreenStatus,
)
from commec.control_list.containers import ControlListOutput, ListMode
from commec.utils.logger import setup_console_logging



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
        "control_hash": 11234,
    }
    row.update(overrides)
    return row


def _make_df(*rows) -> pd.DataFrame:
    """Build a DataFrame from one or more row dicts (using _blast_row defaults)."""
    return pd.DataFrame([_blast_row(**r) for r in rows])


def _mock_control_list(
    status: ListMode = ListMode.COMPLIANCE, acronym: str = "DTL"
) -> MagicMock:
    cl = MagicMock()
    cl.status = status
    cl.acronym = acronym
    cl.display_name = f"{acronym} Display Name"
    return cl


def _mock_control_output(
    name: str = "Dangerous Virus",
    category: str = "Viruses",
    species: str = "Virus sp.",
    genus: str = "Virusgenus",
    list_str: str = "DTL",
    source_text: str = "List entry source",
) -> ControlListOutput:
    out = ControlListOutput()
    out.name = name
    out.category = category
    out.species = species
    out.genus = genus
    out.list = list_str
    out.source_text = source_text
    return out


def _make_screen_result(query_name: str = "seq1") -> ScreenResult:
    result = ScreenResult()
    result.queries[query_name] = QueryResult(query=query_name)
    return result


def _make_match_range() -> MatchRange:
    return MatchRange(1e-10, 100, 300)


# ── _create_hit_info ───────────────────────────────────────────────────────────

class TestCreateHitInfo:

    def _row(self) -> pd.Series:
        return pd.Series(_blast_row())

    def test_basic_fields_populated_without_control_output(self):
        row = self._row()
        result = _create_hit_info(row)
        assert result["evalue"] == row["evalue"]
        assert result["taxid"] == row["subject tax ids"]
        assert result["target_hit"] == row["subject acc."]
        assert result["target_description"] == row["subject title"]
        assert result["query_start"] == row["q. start"]
        assert result["query_end"] == row["q. end"]
        assert result["match_start"] == row["s. start"]
        assert result["match_end"] == row["s. end"]

    def test_no_taxonomy_fields_without_control_output(self):
        result = _create_hit_info(self._row())
        assert "species" not in result
        assert "genus" not in result
        assert "controlled_by_lists" not in result

    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_with_control_output_populates_taxonomy_fields(self, mock_gcl):
        mock_gcl.return_value = _mock_control_list()
        co = _mock_control_output(species="Ebola virus", genus="Ebolavirus")
        result = _create_hit_info(self._row(), [co])
        assert result["species"] == "Ebola virus"
        assert result["genus"] == "Ebolavirus"

    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_with_control_output_populates_controlled_by_lists(self, mock_gcl):
        mock_gcl.return_value = _mock_control_list()
        co = _mock_control_output(source_text="DTL entry for Ebola")
        result = _create_hit_info(self._row(), [co])
        assert len(result["controlled_by_lists"]) == 1
        assert result["controlled_by_lists"][0]["source"] == "DTL entry for Ebola"


# ── _create_hit_result_from_annotations ───────────────────────────────────────

class TestCreateHitResultFromAnnotations:
    """
    Tests the outcome status and hit metadata produced from pre-built annotation
    dicts (the dicts come from _create_hit_info in production).
    """

    def _reg_ann(self, taxid: int = 111, species: str = "Virus sp.") -> dict:
        return {
            "evalue": 1e-10,
            "taxid": taxid,
            "species": species,
            "target_hit": "NC_001",
            "target_description": "Dangerous",
        }

    def _non_reg_ann(self, taxid: int = 999) -> dict:
        return {
            "evalue": 1e-5,
            "taxid": taxid,
            "target_hit": "NC_999",
            "target_description": "Harmless",
        }

    def test_compliance_mode_gives_flag(self):
        co = _mock_control_output()
        hit, _ = _create_hit_result_from_annotations("12345",
            [self._reg_ann()], [], [co], ListMode.COMPLIANCE,
            ScreenStep.TAXONOMY_AA, _make_match_range(),
        )
        assert hit.recommendation.status == ScreenStatus.FLAG

    def test_compliance_warn_mode_gives_warn(self):
        co = _mock_control_output()
        hit, _ = _create_hit_result_from_annotations("12345",
            [self._reg_ann()], [], [co], ListMode.COMPLIANCE_WARN,
            ScreenStep.TAXONOMY_AA, _make_match_range(),
        )
        assert hit.recommendation.status == ScreenStatus.WARN

    def test_conditional_num_mode_gives_pass(self):
        co = _mock_control_output()
        hit, _ = _create_hit_result_from_annotations("12345",
            [self._reg_ann()], [], [co], ListMode.CONDITIONAL_NUM,
            ScreenStep.TAXONOMY_AA, _make_match_range(),
        )
        assert hit.recommendation.status == ScreenStatus.PASS

    def test_non_empty_non_regulated_overrides_status_to_pass(self):
        co = _mock_control_output()
        hit, _ = _create_hit_result_from_annotations("12345",
            [self._reg_ann()], [self._non_reg_ann()], [co], ListMode.COMPLIANCE,
            ScreenStep.TAXONOMY_AA, _make_match_range(),
        )
        assert hit.recommendation.status == ScreenStatus.PASS

    def test_hit_name_uses_longest_compliance_display_name(self):
        co_short = _mock_control_output(name="Short")
        co_long = _mock_control_output(name="A Much Longer Name")
        hit, _ = _create_hit_result_from_annotations("12345",
            [self._reg_ann()], [], [co_short, co_long], ListMode.COMPLIANCE,
            ScreenStep.TAXONOMY_AA, _make_match_range(),
        )
        assert hit.name == "A Much Longer Name"

    def test_log_message_is_non_empty_string(self):
        co = _mock_control_output()
        _, log = _create_hit_result_from_annotations("12345",
            [self._reg_ann()], [], [co], ListMode.COMPLIANCE,
            ScreenStep.TAXONOMY_AA, _make_match_range(),
        )
        assert isinstance(log, str) and len(log) > 0

    def test_step_set_on_hit_recommendation(self):
        co = _mock_control_output()
        hit, _ = _create_hit_result_from_annotations("12345",
            [self._reg_ann()], [], [co], ListMode.COMPLIANCE,
            ScreenStep.TAXONOMY_NT, _make_match_range(),
        )
        assert hit.recommendation.from_step == ScreenStep.TAXONOMY_NT


# ── _create_hit_result_for_cluster ────────────────────────────────────────────

class TestCreateHitResultForCluster:
    """
    Tests the cluster-level hit creation: control-list mode determines flag/warn/pass/none.
    Patches get_regulation (returns ControlListOutput list) and get_control_lists
    (returns a ControlList mock with a .status ListMode).
    """

    def _cluster(self, taxid: int = 111) -> pd.DataFrame:
        return pd.DataFrame([_blast_row(**{"subject tax ids": taxid, "regulated": True})])

    def _non_reg(self) -> pd.DataFrame:
        return pd.DataFrame([],columns=_blast_row().keys())

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_compliance_mode_returns_flag_hit(self, mock_gcl, mock_gr):
        co = _mock_control_output()
        mock_gr.return_value = [co]
        mock_gcl.return_value = _mock_control_list(ListMode.COMPLIANCE)
        hit, _ = _create_hit_result_for_cluster(
            self._cluster(), pd.DataFrame([]), 100, 300, 111, ScreenStep.TAXONOMY_AA
        )
        assert hit is not None
        assert hit.recommendation.status == ScreenStatus.FLAG, f"Expected a ScreenStatus FLAG: \n {json.dumps(asdict(hit), indent=2)}"

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_compliance_warn_mode_returns_warn_hit(self, mock_gcl, mock_gr):
        co = _mock_control_output()
        mock_gr.return_value = [co]
        mock_gcl.return_value = _mock_control_list(ListMode.COMPLIANCE_WARN)
        hit, _ = _create_hit_result_for_cluster(
            self._cluster(), self._non_reg(), 100, 300, 111, ScreenStep.TAXONOMY_AA
        )
        assert hit is not None
        assert hit.recommendation.status == ScreenStatus.WARN

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_conditional_num_below_threshold_returns_none(self, mock_gcl, mock_gr):
        setup_console_logging(logging.DEBUG)
        # N_NON_REGIONAL_HITS_TO_WARN == 2; one taxid → len(conditional_compliances) == 1 < 2
        co = _mock_control_output()
        mock_gr.return_value = [co]
        mock_gcl.return_value = _mock_control_list(ListMode.CONDITIONAL_NUM)
        hit, _ = _create_hit_result_for_cluster(
            self._cluster(taxid=111), self._non_reg(), 100, 300, 111, ScreenStep.TAXONOMY_AA
        )
        assert hit is None

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_conditional_num_at_threshold_returns_pass_hit(self, mock_gcl, mock_gr):
        # Two distinct taxids → len(conditional_compliances) == 2 == N_NON_REGIONAL_HITS_TO_WARN
        co = _mock_control_output()
        mock_gr.return_value = [co]
        mock_gcl.return_value = _mock_control_list(ListMode.CONDITIONAL_NUM)
        cluster = pd.DataFrame([
            _blast_row(**{"subject tax ids": 111, "regulated": True}),
            _blast_row(**{"subject tax ids": 222, "regulated": True, 
                          "subject acc.": "NC_002", "subject title": "Unique Title"}),
        ])
        hit, _ = _create_hit_result_for_cluster(
            cluster, self._non_reg(), 100, 300, 111, ScreenStep.TAXONOMY_AA
        )
        assert hit is not None
        assert hit.recommendation.status == ScreenStatus.PASS, f"Expected a ScreenStatus PASS: \n {json.dumps(asdict(hit))}"

    @patch("commec.screeners.check_reg_path.get_regulation")
    @patch("commec.screeners.check_reg_path.get_control_lists")
    def test_unrecognised_list_mode_returns_none(self, mock_gcl, mock_gr):
        co = _mock_control_output()
        mock_gr.return_value = [co]
        # Status doesn't match COMPLIANCE / COMPLIANCE_WARN / CONDITIONAL_NUM
        mock_gcl.return_value = MagicMock(status=None)
        hit, _ = _create_hit_result_for_cluster(
            self._cluster(), self._non_reg(), 100, 300, 111, ScreenStep.TAXONOMY_AA
        )
        assert hit is None


# ── _get_hit_result_from_data ─────────────────────────────────────────────────

class TestGetHitResultFromData:
    """
    Tests per-query aggregation: top-hit filtering, cluster detection, and
    HitResult accumulation. Inner helpers are patched to isolate this layer.
    """

    @patch("commec.screeners.check_reg_path.find_clusters")
    @patch("commec.screeners.check_reg_path.get_controlled_labels")
    @patch("commec.screeners.check_reg_path.get_top_hits")
    def test_no_regulated_hits_returns_empty_list(self, mock_gth, mock_gcl, mock_fc):
        df = _make_df({"regulated": False})
        mock_gth.return_value = df
        mock_gcl.return_value = df
        hits, logs = _get_hit_result_from_data(df, ScreenStep.TAXONOMY_AA)
        assert hits == []
        assert logs == []

    @patch("commec.screeners.check_reg_path._create_hit_result_for_cluster")
    @patch("commec.screeners.check_reg_path.find_clusters")
    @patch("commec.screeners.check_reg_path.get_controlled_labels")
    @patch("commec.screeners.check_reg_path.get_top_hits")
    def test_regulated_hit_propagates_to_output(self, mock_gth, mock_gcl, mock_fc, mock_chrc):
        df = _make_df({"regulated": True, "subject tax ids": 111})
        mock_gth.return_value = df
        mock_gcl.return_value = df
        mock_fc.return_value = (df.assign(cluster=0), [(100, 300)])
        expected_hit = HitResult(
            HitScreenStatus(ScreenStatus.FLAG, ScreenStep.TAXONOMY_AA),
            "Dangerous Virus", "controlled Viruses", _make_match_range(), {},
        )
        mock_chrc.return_value = (expected_hit, "log line")
        hits, logs = _get_hit_result_from_data(df, ScreenStep.TAXONOMY_AA)
        assert len(hits) == 1
        assert hits[0].recommendation.status == ScreenStatus.FLAG

    @patch("commec.screeners.check_reg_path._create_hit_result_for_cluster")
    @patch("commec.screeners.check_reg_path.find_clusters")
    @patch("commec.screeners.check_reg_path.get_controlled_labels")
    @patch("commec.screeners.check_reg_path.get_top_hits")
    def test_none_hit_from_cluster_not_included(self, mock_gth, mock_gcl, mock_fc, mock_chrc):
        df = _make_df({"regulated": True, "subject tax ids": 111})
        mock_gth.return_value = df
        mock_gcl.return_value = df
        mock_fc.return_value = (df.assign(cluster=0), [(100, 300)])
        mock_chrc.return_value = (None, "")
        hits, logs = _get_hit_result_from_data(df, ScreenStep.TAXONOMY_AA)
        assert hits == []


# ── _build_log_message ─────────────────────────────────────────────────────────

class TestBuildLogMessage:

    def test_output_contains_status_and_taxids(self):
        msg = _build_log_message(
            ScreenStatus.FLAG, "Viruses",
            [MatchRange(1e-10, 100, 300)],
            "controlled Viruses",
            ["111", "222"], [], ["Virus sp."],
        )
        assert ScreenStatus.FLAG in msg
        assert "111" in msg
        assert "222" in msg

    def test_non_reg_taxids_appear_when_present(self):
        msg = _build_log_message(
            ScreenStatus.WARN, "Viruses",
            [MatchRange(1e-10, 100, 300)],
            "mixed Viruses",
            ["111"], ["9999"], ["Virus sp."],
        )
        assert "9999" in msg

    def test_non_reg_taxids_absent_when_empty(self):
        msg = _build_log_message(
            ScreenStatus.FLAG, "Viruses",
            [MatchRange(1e-10, 100, 300)],
            "controlled Viruses",
            ["111"], [], ["Virus sp."],
        )
        assert "9999" not in msg


# ── parse_taxonomy_hits (integration) ─────────────────────────────────────────

class TestParseTaxonomyHits:
    """
    Integration tests for the main entry point. Filesystem, NCBI taxonomy, and
    control-list calls are patched; the BLAST pipeline is exercised through a
    fabricated DataFrame.
    """

    STEP = ScreenStep.TAXONOMY_AA

    def _make_handler(self, has_hits: bool = True) -> MagicMock:
        h = MagicMock()
        h.validate_output.return_value = True
        h.has_hits.return_value = has_hits
        h.out_file = "/mock/out.blast"
        return h

    def test_returns_1_on_invalid_inputs(self):
        handler = MagicMock()
        handler.validate_output.return_value = False
        handler.out_file = "/mock/out"
        result = ScreenResult()
        with patch("os.path.exists", return_value=True):
            ret = parse_taxonomy_hits(
                handler, result, {}, self.STEP, 1,
            )
        assert ret == 1

    def test_returns_0_and_pass_when_no_hits(self):
        handler = self._make_handler(has_hits=False)
        data = _make_screen_result("seq1")
        with patch("os.path.exists", return_value=True):
            ret = parse_taxonomy_hits(
                handler, data, {}, self.STEP, 1,
            )
        assert ret == 0
        assert data.queries["seq1"].status.protein_taxonomy == ScreenStatus.PASS

    def test_returns_0_and_pass_when_no_regulated_hits(self):
        handler = self._make_handler(has_hits=True)
        data = _make_screen_result("seq1")
        blast_df = _make_df({"query acc.": "seq1", "subject tax ids": 99999, "control_hash" : 12345, "regulated" : False})
        with (
            patch("os.path.exists", return_value=True),
            patch("commec.screeners.check_reg_path.read_blast", return_value=blast_df),
            patch("commec.screeners.check_reg_path.is_regulated", return_value=False),
            patch("commec.screeners.check_reg_path.get_controlled_labels", return_value=blast_df),
        ):
            ret = parse_taxonomy_hits(
                handler, data, {}, self.STEP, 1,
            )
        assert ret == 0
        assert data.queries["seq1"].status.protein_taxonomy == ScreenStatus.PASS

    def test_flags_query_with_regulated_hit(self):
        setup_console_logging(logging.DEBUG)
        handler = self._make_handler(has_hits=True)
        data = _make_screen_result("seq1")
        blast_df = _make_df({"query acc.": "seq1", "subject tax ids": 12345, "control_hash" : 12345})
        flag_hit = HitResult(
            HitScreenStatus(ScreenStatus.FLAG, self.STEP),
            "Dangerous Virus", "controlled Viruses Dangerous Virus",
            _make_match_range(), {},
        )
        with (
            patch("os.path.exists", return_value=True),
            patch("commec.screeners.check_reg_path.read_blast", return_value=blast_df),
            patch("commec.screeners.check_reg_path.is_regulated", return_value=True),
            patch("commec.screeners.check_reg_path.get_controlled_labels", return_value=blast_df),
            patch(
                "commec.screeners.check_reg_path._get_hit_result_from_data",
                return_value=([flag_hit], ["log"]),
            ),
        ):
            ret = parse_taxonomy_hits(
                handler, data, {"seq1": MagicMock()}, self.STEP, 1,
            )
        assert ret == 0
        print(f"Expecting a Flag result for protein taxonomy: {json.dumps(asdict(data),indent=2)}")
        assert data.queries["seq1"].hits[0].recommendation.status == ScreenStatus.FLAG, f"Expected a Flag result for protein taxonomy: {
            json.dumps(asdict(data),indent=2)}"

    def test_unknown_query_acc_in_blast_is_skipped_gracefully(self):
        """blast contains a query accession not present in the ScreenResult queries."""
        handler = self._make_handler(has_hits=True)
        data = _make_screen_result("seq1")
        blast_df = _make_df({"query acc.": "unknown_seq", "subject tax ids": 111})
        with (
            patch("os.path.exists", return_value=True),
            patch("commec.screeners.check_reg_path.read_blast", return_value=blast_df),
            patch("commec.screeners.check_reg_path.is_regulated", return_value=True),
            patch(
                "commec.screeners.check_reg_path._get_hit_result_from_data",
                return_value=([], []),
            ),
        ):
            ret = parse_taxonomy_hits(
                handler, data, {}, self.STEP, 1,
            )
        assert ret == 0
        # seq1 was set to PASS by default; nothing should have changed it
        assert data.queries["seq1"].status.protein_taxonomy == ScreenStatus.PASS
