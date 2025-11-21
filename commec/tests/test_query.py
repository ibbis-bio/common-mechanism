from io import StringIO
import os
import pandas as pd
import pytest
import textwrap
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from commec.config.query import Query, QueryTranslation

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")

def test_get_frame_length():
    # 11 nt query
    query = Query(SeqRecord(Seq("atgtgccatgg"), id="test"))
    assert 9 == query._get_frame_length(frame_offset=0)
    assert 9 == query._get_frame_length(frame_offset=1)
    assert 9 == query._get_frame_length(frame_offset=2)

    # 15 nt query
    query = Query(SeqRecord(Seq("atgtgccatggatgc"), id="test"))
    assert 15 == query._get_frame_length(frame_offset=0)
    assert 12 == query._get_frame_length(frame_offset=1)
    assert 12 == query._get_frame_length(frame_offset=2)

    # 16 nt query
    query = Query(SeqRecord(Seq("atgtgccatggatgca"), id="test"))
    assert 15 == query._get_frame_length(frame_offset=0)
    assert 15 == query._get_frame_length(frame_offset=1)
    assert 12 == query._get_frame_length(frame_offset=2)

def test_translate_to_file(tmp_path):
    query = Query(SeqRecord(Seq("atgtgccatgg"), id="test"))

    expected_output = textwrap.dedent(
        """\
        >test_1
        MCH
        >test_2
        CAM
        >test_3
        VPW
        >test_4
        MAH
        >test_5
        HGT
        >test_6
        PWH
        """
    )

    aa_output = tmp_path / "test_translated.faa"

    query.translate(aa_output)

    # Check if the output file exists
    assert aa_output.exists()

    actual_output = aa_output.read_text()
    assert expected_output.strip() == actual_output.strip()


def test_translate():
    """
    Test translation from nucleotide to 6 frames of protein sequences.
    """
    # 11nt query
    query = Query(SeqRecord(Seq("atgtgccatgg"), id="test"))

    # Input sequence: atgtgccatgg
    # Translations:
    # Frame   Pos     Codon split         Translation
    # 1       0       atg tgc cat gg      MCH
    # 2       1       a tgt gcc atg g     CAM
    # 3       2       at gtg cca tgg      VPW
    # 4       -0      cc atg gca cat      MAH
    # 5       -1      c cat ggc aca t     HGT
    # 6       -2      cca tgg cac at      PWH
    expected_translations = [
        QueryTranslation(frame=1, sequence="MCH"),
        QueryTranslation(frame=2, sequence="CAM"),
        QueryTranslation(frame=3, sequence="VPW"),
        QueryTranslation(frame=4, sequence="MAH"),
        QueryTranslation(frame=5, sequence="HGT"),
        QueryTranslation(frame=6, sequence="PWH"),
    ]

    query._translate()
    assert expected_translations == query.translations

    # 15nt query
    query = Query(SeqRecord(Seq("acgcacctgatcgct"), id="test"))


    # Input sequence: acgcacctgatcgct
    # Translations:
    # Frame   Pos     Codon split              Translation
    # 1       0       acg cac ctg atc gct      THLIA
    # 2       1       a cgc acc tga tcg ct     RTXS
    # 3       2       ac gca cct gat cgc t      APDR
    # 4       -0      agc gat cag gtg cgt      SDQVR
    # 5       -1      a gcg atc agg tgc gt     RSGA
    # 6       -2      ag cga tca ggt gcg t     AIRC
    expected_translations = [
        QueryTranslation(frame=1, sequence="THLIA"),
        QueryTranslation(frame=2, sequence="RTXS"),
        QueryTranslation(frame=3, sequence="APDR"),
        QueryTranslation(frame=4, sequence="SDQVR"),
        QueryTranslation(frame=5, sequence="RSGA"),
        QueryTranslation(frame=6, sequence="AIRC"),
    ]

    query._translate()
    assert expected_translations == query.translations


def test_ambigious():
    """
    Test translation from nucleotide to 6 frames of protein sequences using ambigious nts
    | --------------------------------------------------------------------- |
    | Code             | Bases Represented | Meaning                        |
    | ---------------- | ----------------- | ------------------------------ |
    | **A**            | A                 | Adenine                        |
    | **C**            | C                 | Cytosine                       |
    | **G**            | G                 | Guanine                        |
    | **T** (or **U**) | T (or U in RNA)   | Thymine (or Uracil)            |
    | **R**            | A or G            | puRine                         |
    | **Y**            | C or T            | pYrimidine                     |
    | **S**            | G or C            | Strong interaction (3 H-bonds) |
    | **W**            | A or T            | Weak interaction (2 H-bonds)   |
    | **K**            | G or T            | Keto                           |
    | **M**            | A or C            | aMino                          |
    | **B**            | C or G or T       | not A                          |
    | **D**            | A or G or T       | not C                          |
    | **H**            | A or C or T       | not G                          |
    | **V**            | A or C or G       | not T                          |
    | **N**            | A or C or G or T  | any base (completely unknown)  |
    | --------------------------------------------------------------------- |
    | Codon | Ambiguity | Expansions         | Amino Acid |
    | ----- | --------- | ------------------ | ---------- |
    | AAR   | R=A/G     | AAA, AAG           | Lys        |
    | TAY   | Y=C/T     | TAC, TAT           | Tyr        |
    | GCN   | N=A/C/G/T | GCA, GCC, GCG, GCT | Ala        |
    | AAY   | Y=C/T     | AAC, AAT           | Asn        |
    | GAR   | R=A/G     | GAA, GAG           | Glu        |
    | ACM   | M=A/C     | ACA, ACC           | Thr        |
    | CCS   | S=C/G     | CCC, CCG           | Pro        |
    | GGW   | W=A/T     | GGA, GGT           | Gly        |
    | GTK   | K=G/T     | GTG, GTT           | Val        |
    | ATH   | H=A/C/T   | ATA, ATC, ATT      | Ile        |
    | GTD   | D=A/G/T   | GTA, GTG, GTT      | Val        |
    | CCB   | B=C/G/T   | CCC, CCG, CCT      | Pro        |
    | GTV   | V=A/C/G   | GTA, GTC, GTG      | Val        |
    """
    # 11nt query
    #query = Query(SeqRecord(Seq("atntnccatgg"), id="test"))
    #query = Query(SeqRecord(Seq("ATGAARTAYGCNAAYGARACNABGGADCAHGAVACNTGG"), id="test"))
    query = Query(SeqRecord(Seq("ATGAARTAYGCNAAYGARACMCCSGGWGTKATHGTDCCBGTV"), id="test"))

    expected_translations = [
        QueryTranslation(frame=1, sequence="MKYANETPGVIVPV"),
        QueryTranslation(frame=2, sequence="XXXXXXXXXXXXX"),
        QueryTranslation(frame=3, sequence="EXXXXBXXXXXXX"),
        QueryTranslation(frame=4, sequence="XXXXXXXXXXXXXH"),
        QueryTranslation(frame=5, sequence="XXXXXXXXXXXXS"),
        QueryTranslation(frame=6, sequence="TGTITPGVSXAYF"),
    ]

    query._translate()
    assert expected_translations == query.translations, query.translations[2:4]


@pytest.mark.parametrize(
    "input, expected_output",
    [
        ("short_query_name","short_query_name"),
        ("longquerynamewithnotokenstospliton","longquerynamewithnotokens"),
        ("long_query_name_with_common_tokens_to_split_on", "long_query_name_with")
    ],
)
def test_query_id_creation(input, expected_output):
    """ Small tests to ensure the ID creation works to create IDS < 25 characters."""
    output = Query.create_id(input)
    assert output == expected_output
