import os
import argparse
import textwrap

from commec.flag import add_args, run

SCREEN_DIR = os.path.join(os.path.dirname(__file__), "test_data")

def test_flag(tmp_path):
    """We are lazily writing tests for a full run of flag instead of unit tests."""
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args([SCREEN_DIR, "-o", str(tmp_path), "-r"])
    run(args)

    # Check if the output file exists
    status_output = tmp_path / "screen_pipeline_status.csv"
    assert status_output.exists()

    expected_status = textwrap.dedent(
        f"""\
        name,filepath,flag,biorisk,protein,nucleotide,low_concern,virus_flag,bacteria_flag,eukaryote_flag,low_concern_protein,low_concern_rna,low_concern_dna,rationale
        FLAG_TEST_01,{SCREEN_DIR}/flag_tests.json,Flag,Flag,Pass,Pass,Flag,False,False,False,False,False,False,-
        FLAG_TEST_02,{SCREEN_DIR}/flag_tests.json,Flag,Pass,Flag,Pass,Flag,True,False,False,False,False,False,-
        FLAG_TEST_03,{SCREEN_DIR}/flag_tests.json,Flag,Pass,Flag,Pass,Flag,False,True,False,False,False,False,-
        FLAG_TEST_04,{SCREEN_DIR}/flag_tests.json,Flag,Pass,Flag,Pass,Flag,False,False,True,False,False,False,-
        FLAG_TEST_05,{SCREEN_DIR}/flag_tests.json,Flag,Pass,Pass,Flag,Flag,True,True,True,False,False,False,-
        FLAG_TEST_06,{SCREEN_DIR}/flag_tests.json,Pass,Pass,Mixed,Pass,Pass,True,False,False,False,False,False,-
        FCTEST1,{SCREEN_DIR}/functional.json,Flag,Flag,Flag,Flag,Flag,True,False,False,True,True,True,"Matches sequence with pathogenic or toxin function, and protein and nucleotide sequence with regulated organisms; as well as virulence factor; as well as flags cleared as common or non-hazardous."
        """
    )
    actual_status = status_output.read_text()
    assert expected_status.strip() == actual_status.strip()