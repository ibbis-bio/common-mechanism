""" 
Unit test for ensuring that the databases are being called without errors.
Will fail if databases have not been installed as expected, with correct versions.
"""
import os
import pytest
from commec.tools.diamond import DiamondHandler
from commec.tools.blastn import BlastNHandler
from commec.tools.blastx import BlastXHandler
from commec.tools.hmmer import HmmerHandler
from commec.tools.cmscan import CmscanHandler
from commec.tools.search_handler import DatabaseValidationError

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs")

databases_to_implement = [
    [DiamondHandler, "nr_dmnd", "nr"],
    [BlastNHandler, "nt_blast", "nt"],
    [BlastXHandler, "nr_blast", "nr"],
    [HmmerHandler, "benign_db", "benign.hmm"],
    [CmscanHandler, "benign_db", "benign.cmscan"],
]


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
    assert new_db.check_output()

    version: str = new_db.get_version_information()
    assert version

    if os.path.isfile(output_file):
        os.remove(output_file)


bad_databases = [
    [DiamondHandler, "nr_dmnd", "bad"],
    [BlastNHandler, "nt_blast", "bad"],
    [BlastXHandler, "nr_blast", "bad"],
    [HmmerHandler, "benign_db", "bad.hmm"],
    [CmscanHandler, "benign_db", "bad.cmscan"],
    [DiamondHandler, "bad", "bad"],
    [BlastNHandler, "bad", "bad"],
    [BlastXHandler, "bad", "bad"],
    [HmmerHandler, "bad", "bad.hmm"],
    [CmscanHandler, "bad", "bad.cmscan"],
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


def test_diamond_job_and_threads_calculations():
    """
    Tests a range of threads, and diamond database sizes,
    for automatically calculating the optimum number of runs,
    and threads per run. Such that no CPU time is wasted.
    """
    handler = DiamondHandler(
        "commec/tests/test_dbs/nr_dmnd/nr",
        "commec/tests/test_data/single_record.fasta",
        "output.test",
    )
    for max_threads in range(1, 25):
        for n_database_files in range(3, 9):
            concurrent_runs, threads_per_run = handler.determine_runs_and_threads(
                max_threads, n_database_files
            )
            assert concurrent_runs * threads_per_run == max_threads
