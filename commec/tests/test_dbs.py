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

INPUT_QUERY = os.path.join(os.path.dirname(__file__),"test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__),"test_dbs")

databases_to_implement = [
    [DiamondHandler,   "nr_dmnd",     "nr.dmnd"],
    [BlastNHandler,    "nt_blast",    "nt"],
    [BlastXHandler,    "nr_blast",    "nr"],
    [HmmerHandler,       "benign_db",   "benign.hmm"],
    [CmscanHandler,    "benign_db",   "benign.cmscan"],
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

    new_db = input_db[0](db_file, INPUT_QUERY, output_file)
    new_db.search()
    assert new_db.check_output()

    version : str = new_db.get_version_information()
    assert version
    
    if os.path.isfile(output_file):
        os.remove(output_file)