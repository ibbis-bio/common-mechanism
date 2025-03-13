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
    [CmscanHandler, "benign_db", "benign.cm"],
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

n_jobs = [
    None, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
]

@pytest.mark.parametrize("input_jobs", n_jobs)
def test_diamond_job_and_threads_calculations(input_jobs):
    """
    Tests a range of threads, and diamond database sizes,
    for automatically calculating the optimum number of runs,
    and threads per run. Such that no CPU time is wasted.
    No specific expected outcomes, but we can check general expectations:
     - Never exceed max threads.
     No specific expected outcomes, but we check general expectations (e.g. never exceeding max_threads)
    """
    handler = DiamondHandler(
        "commec/tests/test_dbs/nr_dmnd/nr",
        "commec/tests/test_data/single_record.fasta",
        "output.test",
    )
    handler.jobs = input_jobs

    for max_threads in range(1, 25):
        for n_database_files in range(3, 9):
            concurrent_runs, threads_per_run = handler.determine_runs_and_threads(
                max_threads, n_database_files
            )

            # If input jobs is provided, we should never exceed max threads.
            assert concurrent_runs * threads_per_run <= max_threads

            # If no number of input jobs is provided:
            # We should ALWAYS use all available threads.
            # We may use less than Max Threads if the remainder is 0 for the no. of database files
            if input_jobs is None:
                assert ((concurrent_runs * threads_per_run == max_threads) or
                        (concurrent_runs * threads_per_run % n_database_files == 0)), f"""
                {concurrent_runs} runs with {threads_per_run} threads. Input settings:
                {max_threads} max threads, {n_database_files} dbs, {input_jobs} input jobs no.
                """


@pytest.mark.parametrize(
    "input_jobs, max_threads, n_database_files, expected_runs, expected_threads",
    [
        (None, 20, 6, 2, 10), # jobs capped by db count, using all threads
        (None, 8, 5, 1, 5),   # jobs capped by db count, not using all threads
        (3, 12, 6, 3, 4),     # jobs=3 --> 3 runs with 4 threads each
        (10, 20, 5, 5, 4),    # jobs=10 > db=5, capped to 5 runs with 4 threads each
        (20, 10, 5, 5, 2),    # jobs=20 > threads=10, capped to 5 runs with 2 threads each
        (10, 4, 5, 4, 1),     # jobs=10 > db, threads, cappted to 4 runs with 1 thread each
    ]
)
def test_diamond_job_and_threads_calculations_parametrized(
    input_jobs, max_threads, n_database_files, expected_runs, expected_threads
):
    """
    Specific test cases for Diamond Jobs.
    """
    handler = DiamondHandler(
        "commec/tests/test_dbs/nr_dmnd/nr",
        "commec/tests/test_data/single_record.fasta",
        "output.test",
    )
    handler.jobs = input_jobs
    concurrent_runs, threads_per_run = handler.determine_runs_and_threads(
                max_threads, n_database_files
            )

    assert concurrent_runs == expected_runs, f"""
        {input_jobs} jobs, {max_threads} threads failed
        {concurrent_runs} for expected ({expected_runs}) concurrent runs.
        """
    assert threads_per_run == expected_threads, f"""
        {input_jobs} jobs, {max_threads} threads failed
        {threads_per_run} for expected ({expected_threads}) threads per run.
    """
