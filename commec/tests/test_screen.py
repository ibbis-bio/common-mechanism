
import os
import json
from dataclasses import asdict
import textwrap
from unittest.mock import patch
import pytest
import pandas as pd
from commec.screen import run, ScreenArgumentParser, add_args
from commec.config.json_io import get_screen_data_from_json, encode_screen_data_to_json
from commec.utils.json_html_output import generate_html_from_screen_data
from commec.config.result import ScreenResult

DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")

def cheat_taxonomy_info(
    blast: pd.DataFrame,
    _regulated_taxids: list[str],
    _vaccine_taxids: list[str],
    _db_path: str | os.PathLike,
    _threads: int,
):
    """ 
    Used to override get_taxonomic_labels()
    Everything is a virus, and 
    anything with 12345 taxid is regulated.
    """
    blast["regulated"] = blast["subject tax ids"] == 12345
    blast["superkingdom"] = "Viruses"
    blast["phylum"] = "phylum"
    blast["genus"] = "genus"
    blast["species"] = "species"
    return blast

def print_tmp_path_contents(tmp_path):
    """ 
    Debug print, to see what was created within tmp_path
    """
    print(f"Contents of {tmp_path}:")
    for path in tmp_path.rglob("*"):  # Recursively list all files and directories
        print(path.relative_to(tmp_path), "->", "DIR" if path.is_dir() else "FILE")

@pytest.mark.filterwarnings("ignore:.*scattermapbox.*is deprecated.*:DeprecationWarning")
def test_functional_screen(tmp_path, request):
    """
    Full test, utilising --resume to ignore database runs,
    thus including pre-file parsing to create db outputs, and input fasta file.

    NOTE: There is an issue with how the log file is generated when compared to non-pytest,
    and for some reason it is always empty, and not able to then be printed to screen. This can
    make debugging difficult.
    """

    # The three input queries cover all reading frame scenarios.
    desc_1 = "FCTEST1" # 600mer
    seq_1 = textwrap.dedent(
        """\
        aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
        """
    )
    hmmer_biorisk_to_parse = textwrap.dedent(
        """\
        #                                                         --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
        # tname    accession  tlen qname        accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
        Toxin1     12345      1000 FCTEST1_1    21345       200    0          500  10.1   1   1         0         0 1116.0  10.1     1    30     2   31    40   160  1.00 LargeAreaFlag
        Toxin1     12345      1000 FCTEST1_1    21345       200    0         1000  10.1   1   1         0         0 1116.0  10.1    10    20    11   21    40   160  1.00 SmallImportantFlag
        Toxin1     12345      1000 FCTEST1_1    21345       200    0          100  10.1   1   1         0         0 1116.0  10.1    15    25    16   26    40   160  1.00 SmallUnimportantTRIM
        Toxin3     12345      1000 FCTEST1_5    21345       200    0          500  10.1   1   1         0         0 1116.0  10.1     1    30     2   31    40   160  1.00 ReverseExample
        Toxin2     12345      1000 FCTEST1_1    21345       200    0         1000  10.1   1   1         0         0 1116.0  10.1    36    63    36   63    30   60  1.00 Warning Example"""
    )
    blastnr_to_parse = textwrap.dedent(
        """\
        #query acc.	title	        subject acc.taxid	evalue	bit score	% identity	    q.len	q.start	q.end	s.len	s. start	s. end
        FCTEST1	ShouldntClear	    NR_HIT_FLAG1	12345	    0.0	    BITSCORE	99.999	    600	    320	    380	    100	    1	        100
        FCTEST1	ShouldClearBySynbio	NR_HIT_FLAG2	12345	    0.0	    BITSCORE	99.999	    600	    410	    490	    100	    1	        100
        FCTEST1	ShouldntClear	    NR_HIT_FLAG3	12345	    0.0	    BITSCORE	99.999	    600	    410	    500	    100	    1	        100
        FCTEST1	ShouldClear	        NR_HIT_FLAG4	12345	    0.0	    BITSCORE	99.999	    600	    310	    370	    100	    1	        100
        FCTEST1	ShouldMixedReg	    NR_HIT_MIXED	12345	    0.0	    BITSCORE	99.999	    600	    340	    390	    100	    1	        100
        FCTEST1	ShouldMixednonReg	NR_HIT_MIXED	12346	    0.0	    BITSCORE	99.999	    600	    340	    390	    100	    1	        100
        FCTEST1	ShouldMixedNonReg	NR_HIT_MIXED	12347	    0.0	    BITSCORE	99.999	    600	    340	    390	    100	    1	        100
        """
    )
    blastnt_to_parse = textwrap.dedent(
        """\
        #query acc.	title	subject acc.taxid	evalue	bit score	% identity	    q.len	q.start 	q.end	s.len	s. start	s. end
        FCTEST1	SUBJECT	NT_HIT_FLAG1	    12345	    0.0	    BITSCORE	99.999	    600	    220	    280	     60	    1	        100
        FCTEST1	SUBJECT	NT_HIT_FLAG2	    12345	    0.0	    BITSCORE	99.999	    600	    110	    190	     80	    1	        100
        FCTEST1	SUBJECT	NT_HIT_FLAG3	    12345	    0.0	    BITSCORE	99.999	    600	    110	    200	     90	    1	        100
        FCTEST1	Main	NT_HIT_MIXED	    12345	    0.0	    BITSCORE	99.999	    600	    310	    390	    100	    1	        100
        FCTEST1	NonRegMixedWithMain	NT_HIT_MIXED2	    12346	    0.0	    BITSCORE	99.999	    600	    310	    390	    100	    1	        100
        """
    )

    #PROTEIN
    benign_hmmscan_to_parse = textwrap.dedent(
        """\
        #                                                         --- full sequence --- -------------- this domain -------------           hmm coord   ali coord   env coord
        # tname    accession        tlen qname        accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias    from    to  from    to  from    to  acc description of target
        Benign1    Benign1          1000  FCTEST1_1     21345     60       0.0   1000  10.1   1   1         0         0    1116.0  10.1    67   123    67   123    67   123  1.00  BenignHMMClear
        """
    )

    #RNA, mdl = target, seq = query
    benign_cmscan_to_parse = textwrap.dedent(
        """\
        #target name         accession query name                accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value  inc description of target
        #------------------- --------- ------------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ ---------  --- ---------------------
        BENIGNRNA            12346     FCTEST1	                 Q1         50	    100      200       50      150 STRAND TRUNC PASS   GC    10   1000       0.0  100    BenignCMTestOutput
        """
    )

    #SYNBIO:
    benign_blastnt_to_parse = textwrap.dedent(
        """\
        #query acc.	title	subject acc.taxid	evalue	bit score	% identity	    q.len	q.start	q.end	    s.len	s. start	s. end
        FCTEST1	BENIGNSYNBIO	BENIGN_SB	210	        0.0	    BITSCORE	99.999	    600	    410	    480	    90	    10	    100
        """
    )

    input_fasta_path = tmp_path / "functional.fasta"
    input_fasta_path.write_text(f">{desc_1}\n{seq_1}\n")
    json_output_path = tmp_path / "functional.output.json"
    desired_json_output_path = os.path.join(os.path.dirname(__file__), "test_data/functional.json")
    html_output_path = os.path.join(os.path.dirname(__file__), "test_data/functional")

    os.mkdir(tmp_path / "output_functional")
    os.mkdir(tmp_path / "input_functional")

    # --RESUME FILES::
    # BIORISK FILES
    biorisk_db_output_path = tmp_path / "output_functional/functional.biorisk.hmmscan"
    biorisk_db_output_path.write_text(hmmer_biorisk_to_parse)

    # TAXONOMY NR FILES
    nr_db_output_path = tmp_path / "output_functional/functional.nr.blastx"
    nr_db_output_path.write_text(blastnr_to_parse)

    # TAXONOMY NT FILES:
    nt_db_output_path = tmp_path / "output_functional/functional.nt.blastn"
    nt_db_output_path.write_text(blastnt_to_parse)

    # BENIGN FILES:
    benign_hmm_output_path = tmp_path / "output_functional/functional.low_concern.hmmscan"
    benign_hmm_output_path.write_text(benign_hmmscan_to_parse)
    benign_cmscan_output_path = tmp_path / "output_functional/functional.low_concern.cmscan"
    benign_cmscan_output_path.write_text(benign_cmscan_to_parse)
    benign_nt_output_path = tmp_path / "output_functional/functional.low_concern.blastn"
    benign_nt_output_path.write_text(benign_blastnt_to_parse)

    # We patch taxonomic labels to avoid making a mini-taxonomy database.
    with patch("commec.screeners.check_reg_path.get_taxonomic_labels", new=cheat_taxonomy_info), patch(
        "sys.argv",
        ["test.py", str(input_fasta_path), "-d", str(DATABASE_DIRECTORY), "-o", str(tmp_path), "--resume"],
    ):
        parser = ScreenArgumentParser()
        add_args(parser)
        args = parser.parse_args()
        run(args)

    print_tmp_path_contents(tmp_path) # Debug print to ensure correct files are located where they should be.
    
    # TESTING BEGINS:
    assert os.path.isfile(json_output_path)

    output_result : ScreenResult = get_screen_data_from_json(json_output_path)
    generate_html_from_screen_data(output_result, html_output_path)

    # If we are writing exemplare data, do it in raw, to test the json_io simultaneously.
    gen_examples = request.config.getoption("--gen-examples")
    if gen_examples:
        encode_screen_data_to_json(output_result, desired_json_output_path)
        #with open(desired_json_output_path, "w", encoding="utf-8") as json_file:
        #    json.dump(tesT, json_file, indent=4)

    test_result : ScreenResult = get_screen_data_from_json(desired_json_output_path)

    print(test_result)

    assert os.path.isfile(str(tmp_path / "functional.screen.log"))

    with open(str(tmp_path / "functional.screen.log"), "r", encoding = "utf-8") as file:
        while True:
            line = file.readline()
            if not line:  # Breaks loop when EOF is reached
                break
            print(line)



    # This WILL be different, mainly the time,
    # But also potentially blast versions etc.
    output_result.commec_info = None
    test_result.commec_info = None

    # Pytest increments the filename version, so ignore the input file.
    output_result.query_info.file = "REPLACED"
    test_result.query_info.file = "REPLACED"


    # Convert both original and retrieved data to dictionaries and compare
    assert asdict(output_result) == asdict(test_result), (
        f"Functional test does not match predicted output, fix code,"
        f" or if new output is expected, run with --gen-examples\n"
        f"Test JSON Reference data: \n{asdict(output_result)}\n"
        f"Test JSON output data: \n{asdict(test_result)}"
    )

