"""
Unit test for ensuring that the databases are being called without errors.
Will fail if databases have not been installed as expected, with correct versions.
"""

import os
import pytest
from commec.tools.blastn import BlastNHandler
from commec.tools.blastx import BlastXHandler
from commec.tools.hmmer import HmmerHandler
from commec.tools.cmscan import CmscanHandler
from commec.tools.search_handler import DatabaseValidationError

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs")

databases_to_implement = [
    [BlastNHandler, "best_match/nucleotide/", "core_nt"],
    [BlastXHandler, "best_match/protein/", "nr"],
    [HmmerHandler, "low_concern/protein", "benign.hmm"],
    [CmscanHandler, "low_concern/rna", "benign.cm"],
]


def print_tmp_path_contents(tmp_path):
    print(f"Contents of {tmp_path}:")
    for path in tmp_path.rglob("*"):  # Recursively list all files and directories
        print(path.relative_to(tmp_path), "->", "DIR" if path.is_dir() else "FILE")


@pytest.mark.parametrize("input_db", databases_to_implement)
def test_database_can_run(input_db):
    """
    Opens a database object on a test database, and runs the test query on it.
    Fails if commec environment is not setup correctly, or if the database object
    defaults are invalid etc.

    Something similar to this would be useful to be run
    instead of --help during the conda recipe checks.
    """

    db_dir = os.path.join(DATABASE_DIRECTORY, input_db[1])
    db_file = os.path.join(db_dir, input_db[2])

    output_file = "db.out"

    new_db = input_db[0](db_file, INPUT_QUERY, output_file, force=True)
    new_db.search()
    assert new_db.validate_output()

    version: str = new_db.get_version_information()
    assert version

    if os.path.isfile(output_file):
        os.remove(output_file)


bad_databases = [
    [BlastNHandler, "best_match/nucleotide/", "bad"],
    [BlastXHandler, "best_match/protein/", "bad"],
    [HmmerHandler, "low_concern_db", "bad.hmm"],
    [CmscanHandler, "low_concern_db", "bad.cmscan"],
    [BlastNHandler, "bad", "bad"],
    [BlastXHandler, "bad", "bad"],
    [HmmerHandler, "bad", "bad.hmm"],
    [CmscanHandler, "bad", "bad.cmscan"],
    # Wrong BLAST alias: directory has `core_nt.*` files but user configured prefix `nt`.
    # Previous glob-based validation matched any file containing 'nt' and silently passed,
    # surfacing a cryptic "No alias or index file found" error from BLAST itself.
    [BlastNHandler, "best_match/nucleotide/", "nt"],
    [BlastXHandler, "best_match/protein/", "n"],
]


@pytest.mark.parametrize("input_db", bad_databases)
def test_database_no_file(input_db):
    """
    Simply ensures that the input databases are failing there validation.
    """
    db_dir = os.path.join(DATABASE_DIRECTORY, input_db[1])
    db_file = os.path.join(db_dir, input_db[2])
    output_file = "db.out"

    try:
        input_db[0](db_file, INPUT_QUERY, output_file)
        assert False
    except DatabaseValidationError:
        assert True
