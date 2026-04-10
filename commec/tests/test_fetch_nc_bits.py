from io import StringIO
import os
import pandas as pd
import pytest
import textwrap
from Bio import SeqIO
from unittest.mock import patch

from commec.tools.fetch_nc_bits import (
    _get_ranges_with_no_hits,
    calculate_noncoding_regions_per_query,
)

from commec.config.screen_io import ScreenIO
from commec.config.result import QueryResult, ScreenStatus
from commec.screen import add_args, ScreenArgumentParser
from commec.tools.blastx import BlastXHandler

DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs")

@pytest.mark.parametrize(
    "hits, query_length, nc_ranges",
    [
        # Two protein hits, no noncoding regions > 50bp
        ([(1, 50), (100, 150), (175, 299)], 300, []),
        # One protein hit, < 50bp nocoding regions on the ends
        ([(50, 251)], 300, []),
        # One protein hit, > 50bp nocoding regions on the ends
        ([(51, 250)], 300, [(1, 50), (251, 300)]),
        # Three protein hits, one noncoding region >50bp
        (
            [(1, 40), (140, 265), (300, 349)],
            300,
            [(41, 139)],
        ),
        # Regression: reported by Synplogen. When a lower-identity hit
        # encompasses a higher-identity one, the coding region is the union of
        # both hits, so for hits [(400, 600), (200, 800)] on a 1000 bp query the
        # only non-coding regions are [1, 199] and [801, 1000].
        # Previously _get_ranges_with_no_hits walked hits pairwise in
        # sorted-by-start order and used hit_ranges[-1][1] (= 600, the end of
        # [400, 600]) as the rightmost covered base, incorrectly returning
        # [(1, 199), (601, 1000)] and re-screening [601, 800] at the nucleotide
        # level even though it is covered by the [200, 800] protein hit.
        ([(400, 600), (200, 800)], 1000, [(1, 199), (801, 1000)]),
        # Fully nested: inner hit is entirely inside the outer one.
        ([(100, 200), (120, 180)], 300, [(1, 99), (201, 300)]),
        # Partial overlap on the right.
        ([(100, 250), (200, 400)], 500, [(1, 99), (401, 500)]),
        # Adjacent hits (no gap between them) should be treated as one coding block.
        ([(60, 150), (151, 240)], 300, [(1, 59), (241, 300)]),
    ],
)
def test_get_ranges_with_no_hits(hits, query_length, nc_ranges):
    """
    Test the BLAST hits are successfully converted into noncoding ranges.
    """

    def _create_mock_blast_df_from(hits, query_length):
        data = {
            "q. start": [hit[0] for hit in hits],
            "q. end": [hit[1] for hit in hits],
            "query length": [query_length] * len(hits),
        }
        df = pd.DataFrame(data)
        return df.reset_index(drop=True)  # This adds a numeric index

    blast_df = _create_mock_blast_df_from(hits, query_length)
    assert _get_ranges_with_no_hits(blast_df) == nc_ranges


def test_fetch_nocoding_regions(tmp_path):
    """Full test, including file parsing."""

    desc_1 = "NC_TEST01"
    desc_2 = "NC_TEST02"
    seq_1 = textwrap.dedent(
        """\
        ggtagttccctaaacttatcattaagcgatcttcatcgtcaggtatctcgattggtgcagcaagagagcggtgattgt
        accgggaaattaagaggtaacgttgctgccaataaagaaactacctttcaaggtttgaccatagccagtggagccaga
        gagtcagaaaaagtatttgctcaaactgtactaagccacgtagcaaatgttgttctaactcaagaagataccgctaag
        ctattgcaaagtacggtaaagcataatttgaataattatgacttaagaagtgtcggcaatggtaat
        """
    )
    seq_2 = textwrap.dedent(
        """\
        atggcacaagtcattaataccaacagcctctcgctgatcactcaaaataatatcaacaagaaccagtctgcgctgtcg
        agttctatcgagcgtctgtcttctggcttgcgtattaacagcgcgaaggatgacgcagcgggtcaggcgattgctaac
        cgtttcacctctaacattaaaggcctgactcaggcggcccgtaacgccaacgacggtatctccgttgcgcagaccacc
        gaaggcgcgctgtccgaaatcaacaacaacttacagcgtgtgcgtgaactgacggtacaggccact
        """
    )

    blast_to_parse = textwrap.dedent(
        """\
        # BLASTX 2.15.0+
        # Query: NC_TEST
        # Database: /root/commec-dbs/mock
        #query acc.	subject title	subject acc.	subject tax ids	evalue	bit score	% identity	query length	q. start	q. end	subject length	s. start	s. end
        # 3 hits found
        NC_TEST01	SUBJECT	SUBJECT_ACC	TAXID	0.0	BITSCORE	99.999	300	101	200	500	1	100
        NC_TEST02	SUBJECT	SUBJECT_ACC	TAXID	0.0	BITSCORE	99.999	300	25	80	500	1	100
        NC_TEST02	SUBJECT	SUBJECT_ACC	TAXID	0.0	BITSCORE	99.999	300	100	300	500	1	100
        """
    )

    expected_output = textwrap.dedent(
        """\
        >NC_TEST01 (1-100) (201-300)
        ggtagttccctaaacttatcattaagcgatcttcatcgtcaggtatctcgattggtgcagcaagagagcggtgattgtaccgggaaattaagaggtaacgaaatgttgttctaactcaagaagataccgctaagctattgcaaagtacggtaaagcataatttgaataattatgacttaagaagtgtcggcaatggtaat
        """
    )

    # Setup Expected files
    input_fasta = tmp_path / "fetch_nc_input.fasta"
    input_fasta.write_text(f">{desc_1}\n{seq_1}\n>{desc_2}\n{seq_2}\n")
    input_blast = tmp_path / "fetch_nc_input.blastx"
    input_blast.write_text(blast_to_parse)

    # Create Dictionary of queries for funciton input.
    with patch(
        "sys.argv",
        ["test.py", "--skip-tx", str(input_fasta), "-d", str(DATABASE_DIRECTORY), "-o", str(tmp_path)],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        screen_io = ScreenIO(parser.parse_args())
        screen_io.setup()

    queries = screen_io.parse_input_fasta()
    for query in queries.values():
        query.result = QueryResult()

    # Setup result handler for function input.
    db_file = os.path.join(DATABASE_DIRECTORY, "nr_blast/nr")
    handler = BlastXHandler(db_file, input_fasta, input_blast, force=True)

    calculate_noncoding_regions_per_query(handler, queries)

    # Generate the non-coding fasta text.
    actual_output = ""
    for query in queries.values():
        if query.result.status.nucleotide_taxonomy == ScreenStatus.SKIP:
            continue
        actual_output += query.get_non_coding_regions_as_fasta()

    assert actual_output.strip() == expected_output.strip()
