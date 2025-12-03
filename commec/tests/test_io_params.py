import pytest
from unittest.mock import patch
import os
import yaml

from commec.config.screen_io import ScreenIO
from commec.cli import ScreenArgumentParser
from commec.screen import add_args
from commec.utils.file_utils import expand_and_normalize

INPUT_QUERY = os.path.join(os.path.dirname(__file__), "test_data/single_record.fasta")
DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")

@pytest.fixture
def expected_defaults():
    return {
        "base_paths": {
            "default": "commec-dbs/"
        },
        "databases": {
            "low_concern": {
                "rna": {"path": "commec-dbs/low_concern/rna/benign.cm"},
                "dna": {"path": "commec-dbs/low_concern/dna/benign.fasta"},
                "protein": {"path": "commec-dbs/low_concern/protein/benign.hmm"},
                "annotations": 'commec-dbs/low_concern/low_concern_annotations.tsv',
                "taxids": "commec-dbs/low_concern/vax_taxids.txt"
            },
            "biorisk": {
                "path": "commec-dbs/biorisk/biorisk.hmm",
                "annotations": 'commec-dbs/biorisk/biorisk_annotations.csv',
                "taxids": "commec-dbs/biorisk/reg_taxids.txt",
            },
            "regulated_nt": {
                "path": "commec-dbs/nt_blast/core_nt"
            },
            "regulated_protein": {
                "blast": {"path": "commec-dbs/nr_blast/nr"},
                "diamond": {"path": "commec-dbs/nr_dmnd/nr.dmnd"}
            },
            "taxonomy": {
                "path": "commec-dbs/taxonomy/",
            },
            "regulated_lists": {
                "path": "commec-dbs/regulated_lists",
                "regions": "all"
            }
        },
        "threads": 1,
        "protein_search_tool": "blastx",
        "skip_taxonomy_search": False,
        "skip_nt_search": False,
        "do_cleanup": False,
        "diamond_jobs": None,
        "force": False,
        "resume": False,
        "verbose": False
    }

@pytest.fixture
def custom_yaml_config():
    return {
        "databases": {
            "biorisk": {
                "taxids" : "custom_path.txt"
            }
        },
        "skip_taxonomy_search": True,
        "force": True,
        "threads": 8,
    }

@pytest.fixture
def expected_updated_from_custom_yaml():
    return {
        "base_paths": {
            "default": "commec-dbs/"
        },
        "databases": {
            "low_concern": {
                "rna": {"path": "commec-dbs/low_concern/rna/benign.cm"},
                "dna": {"path": "commec-dbs/low_concern/dna/benign.fasta"},
                "protein": {"path": "commec-dbs/low_concern/protein/benign.hmm"},
                "annotations": 'commec-dbs/low_concern/low_concern_annotations.tsv',
                "taxids": "commec-dbs/low_concern/vax_taxids.txt"
            },
            "biorisk": {
                "path": "commec-dbs/biorisk/biorisk.hmm",
                "annotations": 'commec-dbs/biorisk/biorisk_annotations.csv',
                "taxids": "custom_path.txt",
            },
            "regulated_nt": {
                "path": "commec-dbs/nt_blast/core_nt"
            },
            "regulated_protein": {
                "blast": {"path": "commec-dbs/nr_blast/nr"},
                "diamond": {"path": "commec-dbs/nr_dmnd/nr.dmnd"}
            },
            "taxonomy": {
                "path": "commec-dbs/taxonomy/",
            },
            "regulated_lists": {
                "path": "commec-dbs/regulated_lists",
                "regions": "all"
            }
        },
        "threads": 8,
        "protein_search_tool": "blastx",
        "skip_taxonomy_search": True,
        "skip_nt_search": False,
        "do_cleanup": False,
        "diamond_jobs": None,
        "force": True,
        "resume": False,
        "verbose": False
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

def test_user_yaml_override(tmp_path, expected_updated_from_custom_yaml, custom_yaml_config):
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
    assert expected_updated_from_custom_yaml == params.config

def test_cli_override(tmp_path, expected_updated_from_custom_yaml, custom_yaml_config):
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
        "--skip-tx", # skip taxonomy
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
    expected_updated_from_custom_yaml["skip_nt_search"] = True
    expected_updated_from_custom_yaml["do_cleanup"] = True
    db_str_to_override = expected_updated_from_custom_yaml["base_paths"]["default"]

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
    
    expected_defaults = recursive_override(
        expected_updated_from_custom_yaml, db_str_to_override, str(tmp_path) + "/"
    )

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



@pytest.mark.parametrize(
    "base_path, low_concern_path, expected_path",
    [
        # Expected (basepath has terminal separator)
        ("commec-test/", "{default}low_concern/rna/test.cm", "commec-test/low_concern/rna/test.cm"),
        # No separators
        ("commec-test", "{default}low_concern/rna/test.cm", "commec-test/low_concern/rna/test.cm"),
        # Subpath has separator
        ("commec-test", "{default}/low_concern/rna/test.cm", "commec-test//low_concern/rna/test.cm"),
        # Double separators
        ("commec-test/", "{default}/low_concern/rna/test.cm", "commec-test//low_concern/rna/test.cm"),
    ],
)
def test_format_config_paths(tmp_path, base_path, low_concern_path, expected_path):
    config_yaml = {
        "base_paths": {
            "default": base_path
        },
        "databases": {
            "low_concern" : {
                "rna" : {
                    "path": low_concern_path
                }
            }
        }
    }
    user_config_path = tmp_path / "user_config.yaml"
    with open(user_config_path, 'w') as f:
        yaml.dump(config_yaml, f)
    
    parser = ScreenArgumentParser()
    add_args(parser)
    args = parser.parse_args([INPUT_QUERY, "--config", str(user_config_path)])
    params = ScreenIO(args)
    
    assert expected_path == params.config["databases"]["low_concern"]["rna"]["path"]


@pytest.mark.parametrize(
    "input_file, prefix_arg, expected_prefix, is_makedirs_called",
    [
        # No prefix - keeps relative path
        ("dir/file.fasta", None, "dir/file", False),
        # Directory prefix - places in dir
        ("./file.fasta", "dir/output/", "dir/output/file", True),
        # Custom prefix - use that directly
        ("dir/file.fasta", "dir/output", "dir/output", False),
        # User directory prefix - places in expanded dir
        ("dir/file.fasta", "~", "~/file", True),
        # Relative directory prefix - places in dir
        ("dir/file.fasta", "..", "../file", True),
    ],
)
@patch("os.makedirs")
def test_get_output_prefix(
    mock_makedirs, input_file, prefix_arg, expected_prefix, is_makedirs_called
):
    prefix, output_prefix, input_prefix = ScreenIO._get_output_prefixes(input_file, prefix_arg)
    assert expected_prefix == prefix, f"Expected: {expected_prefix}, got {prefix}"

    # Verify makedirs was called when appropriate
    #if is_makedirs_called:
    #    mock_makedirs.assert_called_once_with(expand_and_normalize(prefix_arg), exist_ok=True)
    #else:
    #    mock_makedirs.assert_not_called()