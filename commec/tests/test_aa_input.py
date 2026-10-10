"""
Integration tests for amino acid (protein) query input support in COMMEC.

Covers the protein screening path end-to-end using ScreenTesterFactory:
  - biorisk hmmscan runs directly on protein input (no 6-frame translation)
  - blastp replaces blastx for taxonomy search
  - blastn and RNA low-concern steps are skipped for protein-only queries
"""

from unittest.mock import patch

from commec.config.result import ScreenStatus, ScreenStep
from commec.tests.screen_factory import ScreenTesterFactory


def test_protein_query_no_hits_passes(tmp_path):
    """Protein query with no biorisk or taxonomy hits should report PASS."""
    factory = ScreenTesterFactory("prot_no_hits", tmp_path)
    # Step 6: use is_protein=True so the factory writes protein chars and routes to blastp
    factory.add_query("prot1", 60, is_protein=True)
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.PASS


def test_protein_query_biorisk_flag(tmp_path):
    """Protein query with a regulated biorisk hit should be FLAG."""
    factory = ScreenTesterFactory("prot_biorisk_flag", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    # Step 6: biorisk hits for protein queries use query_name without frame suffix
    factory.add_hit(
        ScreenStep.BIORISK,
        "prot1",
        start=1,
        stop=50,
        title="dangerous_toxin",
        accession="TOX001",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.FLAG


def test_protein_query_biorisk_warn(tmp_path):
    """Protein query with an unregulated biorisk hit should be WARNING."""
    factory = ScreenTesterFactory("prot_biorisk_warn", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    # Step 6: regulated=False means the hit is a virulence factor, not a toxin
    factory.add_hit(
        ScreenStep.BIORISK,
        "prot1",
        start=1,
        stop=50,
        title="virulence_factor",
        accession="VIR001",
        regulated=False,
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.WARN


def test_protein_query_with_numeric_suffix_in_name(tmp_path):
    """Protein query named ORF_1 must not cause a KeyError in append_nt_querylength_info.

    The _base_name helper previously stripped the trailing _1 from ORF_1 because it
    looks like an NT frame suffix, producing ORF which is not a key in the queries dict.
    The fix checks the queries dict directly before attempting suffix stripping.
    """
    factory = ScreenTesterFactory("prot_numeric_suffix", tmp_path)
    factory.add_query("ORF_1", 60, is_protein=True)
    factory.add_hit(
        ScreenStep.BIORISK,
        "ORF_1",
        start=1,
        stop=50,
        title="dangerous_toxin",
        accession="TOX001",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["ORF_1"].status.screen_status == ScreenStatus.FLAG


def test_protein_query_low_concern_hmmer_clears_taxonomy_flag(tmp_path):
    """Low-concern HMMER hit must clear a taxonomy flag on protein input.

    Previously, reset_query_statuses for LOW_CONCERN_DNA and LOW_CONCERN_RNA both
    wrote ScreenStatus.SKIP to the shared low_concern field before parse_low_concern_hits
    ran, so the HMMER low-concern step (which does work on protein input) never applied
    its results and taxonomy flags were never cleared.

    Note: biorisk hits are never clearable by the low-concern step (by design);
    only taxonomy hits are eligible for clearance.
    """
    factory = ScreenTesterFactory("prot_lc_clears_tax", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "prot1",
        start=1,
        stop=55,
        title="regulated pathogen protein",
        accession="REG001",
        taxid=12345,
        species="Dangerous species",
        genus="Dangerous",
        superkingdom="Bacteria",
        regulated=True,
    )
    factory.add_hit(
        ScreenStep.LOW_CONCERN_PROTEIN,
        "prot1",
        start=1,
        stop=55,
        title="housekeeping_protein",
        accession="HOUSE001",
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.CLEARED_FLAG


def test_protein_query_regulated_taxonomy_flags(tmp_path):
    """Protein query whose blastp best match is a regulated pathogen should be FLAG."""
    factory = ScreenTesterFactory("prot_tax_flag", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    # Step 6: TAXONOMY_AA hit for a protein query is written to the blastp output file
    factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "prot1",
        start=1,
        stop=55,
        title="regulated pathogen protein",
        accession="REG001",
        taxid=12345,
        species="Dangerous species",
        genus="Dangerous",
        superkingdom="Bacteria",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.FLAG


def test_protein_only_run_skips_blastx(tmp_path):
    """On a protein-only run the nucleotide FASTA is empty, so blastx must not run.

    Regression test: blastx on an empty query file aborts STEP 2 before the
    blastp block runs. The factory pre-writes an (empty) .nr.blastx and runs
    with --resume, so a BlastXHandler.search call here can only come from the
    unguarded code path.
    """
    factory = ScreenTesterFactory("prot_skip_blastx", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    with patch("commec.tools.blastx.BlastXHandler.search") as blastx_search:
        result = factory.run()
    blastx_search.assert_not_called()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.PASS


def test_protein_query_low_concern_clears_with_underscore_name(tmp_path):
    """A protein query named WP_0001 must match its low-concern hmmscan hits.

    The low-concern lookup strips a trailing _<frame> from hmmscan query names.
    Protein records carry no frame suffix, so WP_0001 was truncated to WP, its
    low-concern table came out empty, and the regulated hit was never cleared.
    Underscored names (e.g. RefSeq WP_/NP_/XP_ accessions) are the common case
    for protein input.
    """
    factory = ScreenTesterFactory("prot_lc_underscore", tmp_path)
    factory.add_query("WP_0001", 60, is_protein=True)
    factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "WP_0001",
        start=1,
        stop=55,
        title="regulated pathogen protein",
        accession="REG001",
        taxid=12345,
        species="Dangerous species",
        genus="Dangerous",
        superkingdom="Bacteria",
        regulated=True,
    )
    factory.add_hit(
        ScreenStep.LOW_CONCERN_PROTEIN,
        "WP_0001",
        start=1,
        stop=55,
        title="housekeeping_protein",
        accession="HOUSE001",
    )
    result = factory.run()
    assert result.queries["WP_0001"].status.screen_status == ScreenStatus.CLEARED_FLAG


def test_protein_query_low_concern_partial_coverage_clears(tmp_path):
    """A 45 aa low-concern overlap must clear a protein hit on a 60 aa query.

    The absolute coverage minimum (50) is defined in nucleotides. Applied
    unchanged to amino acid coordinates it makes clearing about 3x stricter
    than for the same sequence submitted as DNA. For protein queries the
    floor is 50 // 3 = 16 aa; the 80% fractional check is unit-free.
    """
    factory = ScreenTesterFactory("prot_lc_coverage", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "prot1",
        start=1,
        stop=55,
        title="regulated pathogen protein",
        accession="REG001",
        taxid=12345,
        species="Dangerous species",
        genus="Dangerous",
        superkingdom="Bacteria",
        regulated=True,
    )
    factory.add_hit(
        ScreenStep.LOW_CONCERN_PROTEIN,
        "prot1",
        start=5,
        stop=50,
        title="housekeeping_protein",
        accession="HOUSE001",
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.CLEARED_FLAG


def test_nucleotide_frame_name_collision_keeps_develop_behaviour(tmp_path):
    """A nucleotide FASTA with both A and A_1 must keep develop behaviour.

    The hmmscan row for frame 1 of query A is named A_1, which is itself a
    query name here. The hmmer.py fix must restore develop coordinate
    handling: the frame stays 1 and ali 1-17 converts to nt 1-51 (the branch
    previously set frame=0 and skipped conversion). Attribution still follows
    develop: the direct name A_1 wins, so the hit attaches to query A_1. That
    attribution is ambiguous at the name level; resolving it needs input-time
    rejection of colliding names, a follow-up outside this PR.
    """
    factory = ScreenTesterFactory("nt_frame_collision", tmp_path)
    factory.add_query("A", 60)
    factory.add_query("A_1", 60)
    # Frame 1 hit on query A: the factory writes the row as "A_1" with
    # ali coordinates 1-17, which frame 1 converts to nt coordinates 1-51.
    factory.add_hit(
        ScreenStep.BIORISK,
        "A",
        start=4,
        stop=51,
        title="dangerous_toxin",
        accession="TOX001",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["A"].status.screen_status == ScreenStatus.PASS
    assert result.queries["A_1"].status.screen_status == ScreenStatus.FLAG
    hit = result.queries["A_1"].hits[0]
    assert (hit.region.query_start, hit.region.query_end) == (1, 51)


def test_nucleotide_blastx_hit_on_digit_suffixed_query_attributed_correctly(tmp_path):
    """blastx reports bare nucleotide query names, so a hit on insert_1 must
    attach to insert_1, not to insert.

    Regression test: blastx reads the nucleotide FASTA (bare names) and
    translates internally, unlike hmmscan, whose input names carry a frame
    suffix. A get_query that strips a digit suffix on a direct hit would
    misattribute the hit to the base name.
    """
    factory = ScreenTesterFactory("nt_blastx_collision", tmp_path)
    factory.add_query("insert", 60)
    factory.add_query("insert_1", 60)
    factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "insert_1",
        start=1,
        stop=51,
        title="regulated pathogen protein",
        accession="REG001",
        taxid=12345,
        species="Dangerous species",
        genus="Dangerous",
        superkingdom="Bacteria",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["insert_1"].status.screen_status == ScreenStatus.FLAG
    assert result.queries["insert"].status.screen_status == ScreenStatus.PASS
