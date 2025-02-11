from io import StringIO
import os
import pandas as pd
import pytest
import textwrap
from Bio import SeqIO
from commec.tools.fetch_nc_bits import (
    _get_ranges_with_no_hits,
    _write_nc_sequences,
    fetch_noncoding_regions,
)


@pytest.mark.parametrize(
    "hits, nc_ranges",
    [
        # Two protein hits, no noncoding regions > 50bp
        ([(1, 50), (100, 150), (175, 299)], []),
        # One protein hit, < 50bp nocoding regions on the ends
        ([(50, 251)], []),
        # One protein hit, > 50bp nocoding regions on the ends
        ([(51, 250)], [(1, 50), (251, 300)]),
        # Three protein hits, one noncoding region >50bp
        (
            [(1, 40), (140, 265), (300, 349)],
            [(41, 139)],
        ),
    ],
)
def test_get_ranges_with_no_hits(hits, nc_ranges):
    """
    Test the BLAST hits are successfully converted into noncoding ranges.
    """

    def _create_mock_blast_df_from(hits):
        data = {
            "q. start": [hit[0] for hit in hits],
            "q. end": [hit[1] for hit in hits],
            "query length": [300] * len(hits),
        }
        df = pd.DataFrame(data)
        return df.reset_index(drop=True)  # This adds a numeric index

    blast_df = _create_mock_blast_df_from(hits)
    assert _get_ranges_with_no_hits(blast_df) == nc_ranges


def test_write_nc_sequences(tmp_path):
    """
    Test that sequences are successfully written to an output file based on noncoding regions.
    """
    desc = "TEST_FETCH_NC"
    seq = "".join([str(i % 10) for i in range(300)])
    record = list(SeqIO.parse(StringIO(f">{desc}\n{seq}\n"), "fasta"))[0]
    nc_ranges = [[5, 59], [91, 170]]
    outfile = tmp_path / "output.fasta"

    _write_nc_sequences(nc_ranges, record, outfile)

    content = outfile.read_text()
    assert f">{desc} 5-59\n{seq[4:59]}\n" in content
    assert f">{desc} 91-170\n{seq[90:170]}\n" in content


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
        NC_TEST02	SUBJECT	SUBJECT_ACC	TAXID	0.0	BITSCORE	99.999	300	275	300	500	1	100
        """
    )

    expected_output = textwrap.dedent(
        """\
        >NC_TEST01 1-100
        ggtagttccctaaacttatcattaagcgatcttcatcgtcaggtatctcgattggtgcagcaagagagcggtgattgtaccgggaaattaagaggtaacg
        >NC_TEST01 201-300
        aaatgttgttctaactcaagaagataccgctaagctattgcaaagtacggtaaagcataatttgaataattatgacttaagaagtgtcggcaatggtaat
        """
    )

    input_fasta = tmp_path / "fetch_nc_input.fasta"
    input_fasta.write_text(f">{desc_1}\n{seq_1}\n>{desc_2}\n{seq_2}\n")
    input_blast = tmp_path / "fetch_nc_input.blastx"
    input_blast.write_text(blast_to_parse)

    nc_output = tmp_path / "fetch_nc_input.blastx.noncoding.fasta"

    fetch_noncoding_regions(str(input_blast), str(input_fasta))

    # Check if the output file exists
    assert nc_output.exists()

    actual_output = nc_output.read_text()
    assert actual_output.strip() == expected_output.strip()
