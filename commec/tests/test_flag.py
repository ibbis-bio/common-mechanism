import pytest
import os
import argparse
import textwrap

from commec.flag import add_args, run

#SCREEN_DIR = os.path.join(os.path.dirname(__file__), "test_data/screen-files")
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
        name,filepath,flag,biorisk,protein,nucleotide,benign,virus_flag,bacteria_flag,eukaryote_flag,benign_protein,benign_rna,benign_dna
        FCTEST1,{SCREEN_DIR}/functional.json,flag,flag,flag,flag,flag,True,False,False,True,True,True
        biorisk-error-2025-02,{SCREEN_DIR}/screen-files/biorisk-error-2025-02.screen,error,error,-,-,-,,,,,,
        fast-mode-2025-02,{SCREEN_DIR}/screen-files/fast-mode-2025-02.screen,pass,warn,skip,skip,skip,,,,,,
        no-hits-2024-06,{SCREEN_DIR}/screen-files/no-hits-2024-06.screen,pass,pass,-,-,skip,,,,,,
        prot-error-2024-08,{SCREEN_DIR}/screen-files/prot-error-2024-08.screen,error,pass,error,-,-,,,,,,
        prot-hit-not-cleared-2024-06,{SCREEN_DIR}/screen-files/prot-hit-not-cleared-2024-06.screen,flag,pass,flag,skip,not-cleared,True,False,False,False,False,False
        prot-mixed-hit-2024-06,{SCREEN_DIR}/screen-files/prot-mixed-hit-2024-06.screen,pass,pass,mix,skip,skip,False,False,False,,,
        prot-multiple-hits-2024-06,{SCREEN_DIR}/screen-files/prot-multiple-hits-2024-06.screen,flag,warn,flag;mix,pass,not-cleared,False,True,False,True,False,False
        prot-nt-hits-cleared-2024-09,{SCREEN_DIR}/screen-files/prot-nt-hits-cleared-2024-09.screen,pass,pass,flag,flag,cleared,False,True,False,False,True,True
        """
    )
    actual_status = status_output.read_text()
    assert expected_status.strip() == actual_status.strip()

    flags_output = tmp_path / "flags.csv"
    assert flags_output.exists()
    expected_flags = textwrap.dedent(
        f"""\
        filename,biorisk,virulence_factor,regulated_virus,regulated_bacteria,regulated_eukaryote,mixed_regulated_and_non_reg,benign
        {SCREEN_DIR}/functional.json,F,P,F,P,P,P,F
        {SCREEN_DIR}/screen-files/biorisk-error-2025-02.screen,Err,Err,-,-,-,-,-
        {SCREEN_DIR}/screen-files/fast-mode-2025-02.screen,P,F,-,-,-,-,-
        {SCREEN_DIR}/screen-files/no-hits-2024-06.screen,P,P,-,-,-,-,-
        {SCREEN_DIR}/screen-files/prot-error-2024-08.screen,P,P,Err,Err,Err,Err,-
        {SCREEN_DIR}/screen-files/prot-hit-not-cleared-2024-06.screen,P,P,F,P,P,P,F
        {SCREEN_DIR}/screen-files/prot-mixed-hit-2024-06.screen,P,P,P,P,P,F,-
        {SCREEN_DIR}/screen-files/prot-multiple-hits-2024-06.screen,P,F,P,F,P,F,F
        {SCREEN_DIR}/screen-files/prot-nt-hits-cleared-2024-09.screen,P,P,P,F,P,P,P
        """
    )
    assert expected_flags.strip() == flags_output.read_text().strip()

    recomendations_output = tmp_path / "flags_recommended.csv"
    assert recomendations_output.exists()
    expected_recomendations = textwrap.dedent(
        f"""\
        {SCREEN_DIR}/functional.json,F
        {SCREEN_DIR}/screen-files/biorisk-error-2025-02.screen,Err
        {SCREEN_DIR}/screen-files/fast-mode-2025-02.screen,P
        {SCREEN_DIR}/screen-files/no-hits-2024-06.screen,P
        {SCREEN_DIR}/screen-files/prot-error-2024-08.screen,Err
        {SCREEN_DIR}/screen-files/prot-hit-not-cleared-2024-06.screen,F
        {SCREEN_DIR}/screen-files/prot-mixed-hit-2024-06.screen,P
        {SCREEN_DIR}/screen-files/prot-multiple-hits-2024-06.screen,F
        {SCREEN_DIR}/screen-files/prot-nt-hits-cleared-2024-09.screen,P
        """
    )
    assert expected_recomendations.strip() == recomendations_output.read_text().strip()
