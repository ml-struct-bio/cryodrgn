"""Test the config module functionality."""

import os
import tempfile
import warnings
import pytest
from datetime import datetime
from unittest.mock import patch, MagicMock
from cryodrgn import config, utils


def test_load_yaml_config(tmpdir):
    """Test loading a yaml configuration file."""
    config_file = tmpdir.join("test_config.yaml")
    test_config = {"model_args": {"zdim": 8}, "version": "0.3.0"}
    utils.save_yaml(test_config, str(config_file))
    
    loaded_config = config.load(str(config_file))
    assert loaded_config == test_config


def test_load_pkl_config(tmpdir):
    """Test loading a pkl configuration file with deprecation warning."""
    config_file = tmpdir.join("test_config.pkl")
    test_config = {"model_args": {"zdim": 8}, "version": "0.2.0"}
    utils.save_pkl(test_config, str(config_file))
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        loaded_config = config.load(str(config_file))
        
        # Check deprecation warning
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert "deprecated" in str(w[0].message)
        
    assert loaded_config == test_config


def test_load_dict_config():
    """Test loading a configuration from a dictionary."""
    test_config = {"model_args": {"zdim": 8}}
    loaded_config = config.load(test_config)
    assert loaded_config is test_config  # Should return the same object


def test_load_unsupported_extension(tmpdir):
    """Test loading a config file with unsupported extension."""
    config_file = tmpdir.join("test_config.txt")
    config_file.write("some content")
    
    with pytest.raises(RuntimeError, match="Unrecognized config extension"):
        config.load(str(config_file))


def test_save_config_defaults(tmpdir):
    """Test saving a config with default filename and version/time addition."""
    test_config = {"model_args": {"zdim": 8}}
    
    with patch('sys.argv', ['cryodrgn', 'train_vae', 'test.mrcs']):
        with patch('cryodrgn.__version__', '0.3.0'):
            saved_file = config.save(test_config, folder=str(tmpdir))
    
    # Check that default filename was used
    assert saved_file == os.path.join(str(tmpdir), "config.yaml")
    assert os.path.exists(saved_file)
    
    # Load and check that version, time, and cmd were added
    loaded_config = utils.load_yaml(saved_file)
    assert loaded_config["version"] == "0.3.0"
    assert loaded_config["cmd"] == ['cryodrgn', 'train_vae', 'test.mrcs']
    assert "time" in loaded_config
    assert isinstance(loaded_config["time"], datetime)
    assert loaded_config["model_args"]["zdim"] == 8


def test_save_config_custom_filename(tmpdir):
    """Test saving a config with custom filename."""
    test_config = {"model_args": {"zdim": 10}}
    custom_filename = "my_config.yaml"
    
    saved_file = config.save(test_config, filename=custom_filename, folder=str(tmpdir))
    
    assert saved_file == os.path.join(str(tmpdir), custom_filename)
    assert os.path.exists(saved_file)


def test_save_config_preserve_existing_metadata():
    """Test that existing version/time/cmd are preserved."""
    existing_time = datetime(2020, 1, 1, 12, 0, 0)
    test_config = {
        "model_args": {"zdim": 6},
        "version": "0.2.0",
        "time": existing_time,
        "cmd": ["old", "command"]
    }
    
    with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as tmp:
        try:
            saved_file = config.save(test_config, filename=tmp.name)
            loaded_config = utils.load_yaml(saved_file)
            
            # Should preserve existing metadata
            assert loaded_config["version"] == "0.2.0"
            assert loaded_config["time"] == existing_time
            assert loaded_config["cmd"] == ["old", "command"]
        finally:
            os.unlink(tmp.name)


def test_update_config_v1():
    """Test updating a v1 config with new parameters."""
    old_config = {
        "model_args": {
            "pe_type": "linear",
            "zdim": 8
        }
    }
    
    updated_config = config.update_config_v1(old_config)
    
    # Should add feat_sigma as None when pe_type is not gaussian
    assert updated_config["model_args"]["feat_sigma"] is None
    assert updated_config["model_args"]["activation"] == "relu"
    assert "tilt_params" in updated_config["model_args"]


def test_update_config_v1_gaussian_pe():
    """Test update_config_v1 with gaussian positional encoding."""
    config_with_gaussian = {
        "model_args": {
            "pe_type": "gaussian",
            "feat_sigma": 0.5,
            "zdim": 8
        }
    }
    
    # Should raise assertion error when feat_sigma is not present with gaussian pe_type
    with pytest.raises(AssertionError):
        old_config = {
            "model_args": {
                "pe_type": "gaussian",
                "zdim": 8
            }
        }
        config.update_config_v1(old_config)


def test_overwrite_config_deprecated():
    """Test the deprecated _overwrite_config function."""
    mock_args = MagicMock()
    mock_args.norm = [0, 1]
    mock_args.D = 128
    mock_args.l_extent = 0.5
    mock_args.qlayers = 3
    mock_args.qdim = 256
    mock_args.zdim = 10
    mock_args.encode_mode = "resid"
    mock_args.players = 3
    mock_args.pdim = 256
    mock_args.enc_mask = None
    mock_args.pe_type = "linear"
    mock_args.feat_sigma = None
    mock_args.pe_dim = None
    mock_args.domain = "fourier"
    mock_args.activation = "relu"
    
    old_config = {
        "dataset_args": {"norm": None},
        "lattice_args": {"D": 65, "extent": 0.25},
        "model_args": {"zdim": 8, "activation": "relu"}
    }
    
    updated_config = config._overwrite_config(old_config, mock_args)
    
    # Check that values were overwritten
    assert updated_config["dataset_args"]["norm"] == [0, 1]
    assert updated_config["lattice_args"]["D"] == 129  # D + 1
    assert updated_config["lattice_args"]["extent"] == 0.5
    assert updated_config["model_args"]["zdim"] == 10
    assert updated_config["model_args"]["qlayers"] == 3


def test_overwrite_config_none_values():
    """Test _overwrite_config with None values (should not overwrite)."""
    mock_args = MagicMock()
    mock_args.norm = None
    mock_args.D = None
    mock_args.zdim = None
    
    old_config = {
        "dataset_args": {"norm": [0, 1]},
        "lattice_args": {"D": 65},
        "model_args": {"zdim": 8}
    }
    
    updated_config = config._overwrite_config(old_config, mock_args)
    
    # Should not overwrite when args are None
    assert updated_config["dataset_args"]["norm"] == [0, 1]
    assert updated_config["lattice_args"]["D"] == 65
    assert updated_config["model_args"]["zdim"] == 8