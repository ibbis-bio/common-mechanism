from io import StringIO
import os
import pandas as pd
import pytest
import textwrap
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from commec.config.query import Query, QueryTranslation

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")


def test_translate(tmp_path):
    """Full test, including file parsing."""
    # Set up query
    query = Query(SeqRecord(Seq("atgtgccatgg"), id="test"))

    expected_translations = [
        QueryTranslation(frame=1, sequence="MCH", nt_start=0, nt_end=9),
        QueryTranslation(frame=2, sequence="CAM", nt_start=1, nt_end=10),
        QueryTranslation(frame=3, sequence="VPW", nt_start=2, nt_end=11),
        QueryTranslation(frame=4, sequence="MAH", nt_start=0, nt_end=9),
        QueryTranslation(frame=5, sequence="HGT", nt_start=1, nt_end=10),
        QueryTranslation(frame=6, sequence="PWH", nt_start=2, nt_end=11),
    ]

    query.translate()
    assert expected_translations == query.translations
