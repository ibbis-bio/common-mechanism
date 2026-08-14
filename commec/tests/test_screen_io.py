import os
from unittest.mock import patch

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from commec.config.screen_io import IoValidationError, ScreenIO, substitute_non_iupac
from commec.screen import ScreenArgumentParser, add_args


@pytest.fixture
def test_data_dir():
    return os.path.join(os.path.dirname(__file__), "test_data")


@pytest.fixture
def database_dir():
    return os.path.join(os.path.dirname(__file__), "test_dbs")


@pytest.mark.parametrize(
    "fasta_name",
    [
        "single_record.fasta",
        "multiple_records.fasta",
        "has_empty_record.fasta",
        "has_empty_description.fasta",
        "has_records_with_same_description.fasta",
    ],
)
def test_default_parameters(fasta_name, test_data_dir, database_dir, tmp_path):
    input_fasta = os.path.join(test_data_dir, fasta_name)
    with patch(
        "sys.argv",
        ["test.py", "--skip-tx", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIO(parser.parse_args())
        assert screen_io.setup()


@pytest.mark.parametrize(
    "fasta_name,expected_record_count",
    [
        pytest.param("single_record.fasta", 1),
        pytest.param("multiple_records.fasta", 2),
    ],
)
def test_parse_input_fasta(
    fasta_name, expected_record_count, test_data_dir, database_dir, tmp_path
):
    input_fasta = os.path.join(test_data_dir, fasta_name)
    with patch(
        "sys.argv",
        ["test.py", "--skip-tx", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIO(parser.parse_args())
        screen_io.setup()

    queries = screen_io.parse_input_fasta()
    assert len(queries) == expected_record_count


@pytest.mark.parametrize(
    "fasta_name",
    [
        "has_empty_record.fasta",
        "has_empty_description.fasta",
        "has_records_with_same_description.fasta",
    ],
)
def test_parse_invalid_input_fasta(fasta_name, test_data_dir, database_dir, tmp_path):
    input_fasta = os.path.join(test_data_dir, fasta_name)
    with patch(
        "sys.argv",
        ["test.py", "--skip-tx", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIO(parser.parse_args())
        screen_io.setup()

    with pytest.raises(IoValidationError):
        screen_io.parse_input_fasta()


@pytest.mark.parametrize(
    "sequence,expected,expected_substitutions",
    [
        pytest.param("ATGCRYWSMKHBVDN", "ATGCRYWSMKHBVDN", 0, id="iupac_codes_kept"),
        # Uracil is a valid base; masking it would destroy RNA queries, and blastn
        # maps 'U' to 'T' natively, as does Biopython when translating
        pytest.param("AUGCUUAAA", "AUGCUUAAA", 0, id="rna_uracil_kept"),
        pytest.param("atgcry", "ATGCRY", 0, id="lower_case_upper_cased_not_masked"),
        pytest.param("atgI*gc?", "ATGNNGCN", 3, id="non_iupac_substituted"),
        # Underscore written by `_write_clean_fasta` in place of a non-ASCII character
        pytest.param("ATG_GCC", "ATGNGCC", 1, id="non_ascii_placeholder_substituted"),
    ],
)
def test_substitute_non_iupac(sequence, expected, expected_substitutions):
    record = SeqRecord(Seq(sequence), id="test_record")
    substitutions = substitute_non_iupac(record)

    assert str(record.seq) == expected
    # Sequence length must be preserved, so that hit coordinates remain valid
    assert len(record.seq) == len(sequence)
    assert substitutions == expected_substitutions


def test_write_clean_fasta_handles_non_ascii_whitespace(
    test_data_dir, database_dir, tmp_path, caplog
):
    input_fasta = os.path.join(test_data_dir, "has_non_ascii_whitespace.fasta")
    with patch(
        "sys.argv",
        ["test.py", "--skip-tx", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIO(parser.parse_args())
        screen_io.setup()

    with open(screen_io.nt_path, "r", encoding="utf-8") as cleaned_fasta:
        header, sequence = cleaned_fasta.read().splitlines()

    # Non-breaking spaces in the header are substituted, so that the record id is preserved
    assert header == ">construct with_non-breaking space"
    # Non-breaking space, vertical tab and form feed are formatting, so they are removed
    assert sequence == (
        "ATGTGCCATGGAATGCGCCATGGAATGCATGGACATGGAATGCATGGAATGCATGGACATGGAATGCATGGA"
    )

    # Whitespace is gone by the time the sequence is parsed, so nothing is masked with 'N'
    with caplog.at_level("WARNING"):
        queries = screen_io.parse_input_fasta()

    assert queries["construct"].sequence == sequence
    assert "substituted" not in caplog.text


def test_write_clean_fasta_handles_non_ascii_characters(
    test_data_dir, database_dir, tmp_path, caplog
):
    input_fasta = os.path.join(test_data_dir, "has_non_ascii_characters.fasta")
    with patch(
        "sys.argv",
        ["test.py", "--skip-tx", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIO(parser.parse_args())
        screen_io.setup()

    with open(screen_io.nt_path, "r", encoding="utf-8") as cleaned_fasta:
        header, sequence = cleaned_fasta.read().splitlines()

    # A byte order mark is consumed while reading, so the header is still a header, and
    # non-breaking, figure and ideographic spaces are all substituted, not just '\xa0'
    assert header == ">construct_with_exotic_whitespace"
    # The en-dash is not whitespace, so it is replaced with an underscore, not removed
    assert "–" not in sequence

    with caplog.at_level("WARNING"):
        queries = screen_io.parse_input_fasta()

    query = queries["construct_with_exotic_whitespace"]
    # Masked with 'N' rather than deleted, so that hit coordinates remain valid
    assert query.sequence == sequence.replace("_", "N")
    assert len(query.sequence) == len(sequence)
    # ...and counted, so that a query full of non-ASCII characters trips the threshold
    assert "substituted 1 non-IUPAC nucleotide characters" in caplog.text


@pytest.mark.parametrize(
    "fasta_name,record_id,expect_not_dna_warning",
    [
        # 2 substitutions in 72 bases is well under NON_IUPAC_SUBSTITUTION_THRESHOLD
        pytest.param(
            "has_non_iupac.fasta",
            "synthetic_construct_with_inosine",
            False,
            id="below_threshold",
        ),
        # Amino acid codes are mostly not IUPAC nucleotide codes, so this blows the threshold
        pytest.param(
            "is_protein_sequence.fasta",
            "haemagglutinin_fragment",
            True,
            id="above_threshold",
        ),
    ],
)
def test_parse_input_fasta_warns_when_substitutions_exceed_threshold(
    fasta_name,
    record_id,
    expect_not_dna_warning,
    test_data_dir,
    database_dir,
    tmp_path,
    caplog,
):
    input_fasta = os.path.join(test_data_dir, fasta_name)
    with patch(
        "sys.argv",
        ["test.py", "--skip-tx", input_fasta, "-d", database_dir, "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIO(parser.parse_args())
        screen_io.setup()

    with caplog.at_level("WARNING"):
        screen_io.parse_input_fasta()

    not_dna_warnings = [
        r.getMessage()
        for r in caplog.records
        if "may not be a nucleotide sequence" in r.getMessage()
    ]

    if not expect_not_dna_warning:
        assert not not_dna_warnings
        return

    assert len(not_dna_warnings) == 1
    assert record_id in not_dna_warnings[0]
