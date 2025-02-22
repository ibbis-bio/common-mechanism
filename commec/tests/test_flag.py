import pytest
import os
import argparse
import textwrap

from commec.flag import add_args, run

SCREEN_DIR = os.path.join(os.path.dirname(__file__), "test_data/screen-files")

def test_flag_status(tmp_path):
    """We are lazily writing tests for a full run of flag instead of unit test."""
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args([SCREEN_DIR, "-o", str(tmp_path)])
    run(args)

    # Check if the output file exists
    status_output = tmp_path / "screen_pipeline_status.csv"
    assert status_output.exists()

    expected_status = textwrap.dedent(
        f"""\
        name,filepath,biorisk,protein,nucleotide,benign
        biorisk-error-2025-02,{SCREEN_DIR}/biorisk-error-2025-02.screen,error,-,-,-
        fast-mode-2025-02,{SCREEN_DIR}/fast-mode-2025-02.screen,warn,skip,skip,skip
        no-hits-2024-06,{SCREEN_DIR}/no-hits-2024-06.screen,pass,-,-,skip
        prot-error-2024-08,{SCREEN_DIR}/prot-error-2024-08.screen,pass,error,-,-
        prot-hit-not-cleared-2024-06,{SCREEN_DIR}/prot-hit-not-cleared-2024-06.screen,pass,flag,skip,not-cleared
        prot-mixed-hit-2024-06,{SCREEN_DIR}/prot-mixed-hit-2024-06.screen,pass,mix,skip,skip
        prot-multiple-hits-2024-06,{SCREEN_DIR}/prot-multiple-hits-2024-06.screen,warn,flag;mix,pass,not-cleared
        prot-nt-hits-cleared-2024-09,{SCREEN_DIR}/prot-nt-hits-cleared-2024-09.screen,pass,flag,flag,cleared
        """
    )
    actual_status = status_output.read_text()
    assert expected_status.strip() == actual_status.strip()

    flags_output = tmp_path / "flags.csv"
    assert flags_output.exists()
    expected_flags = textwrap.dedent(
        f"""\
        filename,biorisk,virulence_factor,regulated_virus,regulated_bacteria,regulated_eukaryote,mixed_regulated_and_non_reg,benign
        {SCREEN_DIR}/biorisk-error-2025-02.screen,Err,Err,-,-,-,-,-
        {SCREEN_DIR}/fast-mode-2025-02.screen,P,F,-,-,-,-,-
        {SCREEN_DIR}/no-hits-2024-06.screen,P,P,-,-,-,-,-
        {SCREEN_DIR}/prot-error-2024-08.screen,P,P,Err,Err,Err,Err,-
        {SCREEN_DIR}/prot-hit-not-cleared-2024-06.screen,P,P,F,P,P,P,F
        {SCREEN_DIR}/prot-mixed-hit-2024-06.screen,P,P,P,P,P,F,-
        {SCREEN_DIR}/prot-multiple-hits-2024-06.screen,P,F,P,F,P,F,F
        {SCREEN_DIR}/prot-nt-hits-cleared-2024-09.screen,P,P,P,F,P,P,P
        """
    )
    assert expected_flags.strip() == flags_output.read_text().strip()

    recomendations_output = tmp_path / "flags_recommended.csv"
    assert recomendations_output.exists()
    expected_recomendations = textwrap.dedent(
        f"""\
        {SCREEN_DIR}/biorisk-error-2025-02.screen,Err
        {SCREEN_DIR}/fast-mode-2025-02.screen,P
        {SCREEN_DIR}/no-hits-2024-06.screen,P
        {SCREEN_DIR}/prot-error-2024-08.screen,Err
        {SCREEN_DIR}/prot-hit-not-cleared-2024-06.screen,F
        {SCREEN_DIR}/prot-mixed-hit-2024-06.screen,P
        {SCREEN_DIR}/prot-multiple-hits-2024-06.screen,F
        {SCREEN_DIR}/prot-nt-hits-cleared-2024-09.screen,P
        """
    )
    assert expected_recomendations.strip() == recomendations_output.read_text().strip()
