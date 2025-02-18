import pytest
from unittest.mock import patch
import os
import yaml

from commec.config.io_parameters import ScreenIO
from commec.cli import ScreenArgumentParser
from commec.screen import add_args

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")

@pytest.fixture
def expected_defaults():
    return {
        "base_paths": {
            "default": "commec-dbs/"
        },
        "databases": {
            "benign": {
                "cm": {"path": "commec-dbs/benign_db/benign.cm"},
                "fasta": {"path": "commec-dbs/benign_db/benign.fasta"},
                "hmm": {"path": "commec-dbs/benign_db/benign.hmm"}
            },
            "biorisk_hmm": {
                "path": "commec-dbs/biorisk_db/biorisk.hmm"
            },
            "regulated_nt": {
                "path": "commec-dbs/nt_blast/nt"
            },
            "regulated_protein": {
                "blast": {"path": "commec-dbs/nr_blast/nr"},
                "diamond": {"path": "commec-dbs/nr_dmnd/nr.dmnd"}
            },
            "taxonomy": {
                "taxonomy_directory": "commec-dbs/taxonomy/",
                "regulated_vaxids": "commec-dbs/biorisk_db/vaxids.txt",
                "benign_taxids": "commec-dbs/benign_db/taxids.txt"
        },
        },
        "threads": 1,
        "protein_search_tool": "blastx",
        "in_fast_mode": False,
        "skip_nt_search": False,
        "do_cleanup": False,
        "diamond_jobs": None,
        "force": False,
        "resume": False
    }

@pytest.fixture
def custom_yaml_config():
    return {
        "databases": {
            "taxonomy": {
                "regulated_vaxids" : "custom_path.txt"
            }
        },
        "in_fast_mode": True,
        "force": True,
        "threads": 8
    }

def test_missing_input_file():
    args = ScreenArgumentParser()
    add_args(args)
    with pytest.raises(SystemExit):
        args = args.parse_args()
    
def test_default_config_only(expected_defaults):
    """Test that default config is loaded when no overrides exist"""
    parser = ScreenArgumentParser()
    add_args(parser)
    args = parser.parse_args([INPUT_QUERY])
    params = ScreenIO(args)
    
    assert expected_defaults == params.config

def test_user_yaml_override(tmp_path, expected_defaults, custom_yaml_config):
    """Test that user YAML properly overrides default config"""
    # Create user config
    user_config_path = tmp_path / "user_config.yaml"
    with open(user_config_path, 'w') as f:
        yaml.dump(custom_yaml_config, f)
    
    parser = ScreenArgumentParser()
    add_args(parser)
    args = parser.parse_args([INPUT_QUERY, "--config", str(user_config_path)])
    params = ScreenIO(args)
    
    # Check that user YAML values override defaults
    expected_defaults.update(custom_yaml_config)

    assert expected_defaults == params.config

def test_cli_override(tmp_path, expected_defaults, custom_yaml_config):
    """Test that CLI args properly override both YAML configs"""
    # Create user config
    user_config_path = tmp_path / "user_config.yaml"
    with open(user_config_path, 'w') as f:
        yaml.dump(custom_yaml_config, f)

    # Add CLI args
    cli_args = [
        INPUT_QUERY,
        "--config",
        str(user_config_path),
        "-f", # fast mode
        "--skip-nt", # skip nt search
        "-c", # do_cleanup
        "-d",
        str(tmp_path)
    ]

    parser = ScreenArgumentParser()
    add_args(parser)
    args = parser.parse_args(cli_args)
    params = ScreenIO(args)
    
    # Override defaults with user YAML
    expected_defaults.update(custom_yaml_config)
    expected_defaults["skip_nt_search"] = True
    expected_defaults["do_cleanup"] = True
    db_str_to_override = expected_defaults["base_paths"]["default"]

    def recursive_override(dictionary, str_to_override, override_str):
        """
        Recursively apply string formatting to read paths from nested yaml config dicts.
        """
        if isinstance(dictionary, dict):
            return {key : recursive_override(value, str_to_override, override_str) 
                    for key, value in dictionary.items()}
        if isinstance(dictionary, str):
           return dictionary.replace(str_to_override, override_str)
        return dictionary
    
    expected_defaults = recursive_override(expected_defaults, db_str_to_override, str(tmp_path))

    assert expected_defaults == params.config

def test_missing_default_config():
    """Test that missing default config raises appropriate error"""
    with patch("importlib.resources.files") as mock_files:
        mock_files.return_value.joinpath.return_value.exists.return_value = False
        args = ScreenArgumentParser()
        add_args(args)
        args = args.parse_args([INPUT_QUERY])
        
        with pytest.raises(FileNotFoundError, match="No default yaml found"):
            _ = ScreenIO(args)

