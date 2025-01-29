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

    query.translate()
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

    query.translate()
    assert expected_translations == query.translations
