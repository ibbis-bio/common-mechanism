from io import StringIO
import os
import pandas as pd
import pytest
import textwrap
from Bio import SeqIO
from commec.config.query import Query

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")


def test_write_six_frame_translation(tmp_path):
    """Full test, including file parsing."""
    translated_faa_path = tmp_path / "test_translate.faa"

    # Set up query
    query = Query(INPUT_QUERY)
    query.nt_path = INPUT_QUERY
    query.aa_path = translated_faa_path
    query.parse_query_data()
    query._write_six_frame_translation()


#    EMBL actual_output
        # >BBa_K380009_1
        # VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPK
        # >BBa_K380009_2
        # XTTNSTKNNKTRSMRSYIYLTXTKNNETPSSKVXKMTQAKALTFXQKLKSXMMLRRRX
        # >BBa_K380009_3
        # RQQIQQRTTKRVLXDLTFTXLKRRTTKRLHPKFKRXPKPKRXPFSRSXKAKXCSGAEX
        # >BBa_K380009_4
        # FRRLSIIXLFSFCXKVSALAWVIFXTLDEGVSLFFVXVRXMXDLIERVLLFFVEFVVY
        # >BBa_K380009_5
        # SAPEHHLAFXLLLKGXRFGLGHLLNFGXRRFVVLRLSXVNVRSHRTRFVVLCXICCLX
        # >BBa_K380009_6
        # FGAXASFSFLASAKRLALWLGSSFKLWMKAFRCSSFKLGKCKISXNAFCCSLLNLLST
#         """

    expected_output = textwrap.dedent(
        """\
        >BBa_K380009_1
        VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPK
        >BBa_K380009_4
        FRRLSIIXLFSFCXKVSALAWVIFXTLDEGVSLFFVXVRXMXDLIERVLLFFVEFVVY
        >BBa_K380009_2
        XTTNSTKNNKTRSMRSYIYLTXTKNNETPSSKVXKMTQAKALTFXQKLKSXMMLRRR
        >BBa_K380009_5
        SAPEHHLAFXLLLKGXRFGLGHLLNFGXRRFVVLRLSXVNVRSHRTRFVVLCXICCL
        >BBa_K380009_3
        RQQIQQRTTKRVLXDLTFTXLKRRTTKRLHPKFKRXPKPKRXPFSRSXKAKXCSGAE
        >BBa_K380009_6
        FGAXASFSFLASAKRLALWLGSSFKLWMKAFRCSSFKLGKCKISXNAFCCSLLNLLST
        """
    )

    # Check if the output file exists
    assert translated_faa_path.exists()

    actual_output = translated_faa_path.read_text()
    assert actual_output.strip() == expected_output.strip()
