"""Test the view_config command functionality."""

import argparse
import os
import tempfile
import warnings
import pytest
from unittest.mock import patch
from cryodrgn.commands import view_config
from cryodrgn import utils


def test_add_args():
    """Test that add_args correctly sets up argument parser."""
    parser = argparse.ArgumentParser()
    view_config.add_args(parser)
    
    # Test that the workdir argument was added
    args = parser.parse_args(["/test/path"])
    assert args.workdir == "/test/path"


def test_main_with_yaml_config(tmpdir):
    """Test main function with a yaml config file."""
    workdir = tmpdir.mkdir("test_workdir")
    config_file = workdir.join("config.yaml")
    
    # Create a test config
    test_config = {
        "version": "0.3.0",
        "time": "2023-01-01 12:00:00",
        "cmd": ["cryodrgn", "train_vae", "test.mrcs"],
        "model_args": {"zdim": 8}
    }
    utils.save_yaml(test_config, str(config_file))
    
    # Mock args
    args = argparse.Namespace(workdir=str(workdir))
    
    # Test that deprecation warning is raised
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with patch('builtins.print') as mock_print:
            view_config.main(args)
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "deprecated" in str(w[0].message)
        
        # Check that the command was printed
        mock_print.assert_called_once_with("cryodrgn train_vae test.mrcs")


def test_main_with_pkl_config(tmpdir):
    """Test main function with a pkl config file."""
    workdir = tmpdir.mkdir("test_workdir")
    config_file = workdir.join("config.pkl")
    
    # Create a test config
    test_config = {
        "version": "0.2.0",
        "model_args": {"zdim": 10}
    }
    utils.save_pkl(test_config, str(config_file))
    
    # Mock args
    args = argparse.Namespace(workdir=str(workdir))
    
    # Test that deprecation warning is raised
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        view_config.main(args)
        assert len(w) >= 1  # At least one warning (could be more due to pkl deprecation)


def test_main_no_config_file(tmpdir):
    """Test main function when no config file exists."""
    workdir = tmpdir.mkdir("test_workdir")
    args = argparse.Namespace(workdir=str(workdir))
    
    with pytest.raises(RuntimeError, match="A configuration file was not found"):
        view_config.main(args)


def test_main_yaml_preferred_over_pkl(tmpdir):
    """Test that yaml config is preferred when both yaml and pkl exist."""
    workdir = tmpdir.mkdir("test_workdir")
    yaml_file = workdir.join("config.yaml")
    pkl_file = workdir.join("config.pkl")
    
    # Create different configs
    yaml_config = {"source": "yaml", "version": "0.3.0"}
    pkl_config = {"source": "pkl", "version": "0.2.0"}
    
    utils.save_yaml(yaml_config, str(yaml_file))
    utils.save_pkl(pkl_config, str(pkl_file))
    
    args = argparse.Namespace(workdir=str(workdir))
    
    # Mock the logger and print to capture output
    with patch('cryodrgn.commands.view_config.logger') as mock_logger:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Ignore deprecation warning
            view_config.main(args)
        
        # Check that yaml version was used
        mock_logger.info.assert_any_call('Version: 0.3.0')


def test_main_minimal_config(tmpdir):
    """Test main function with minimal config (no version, time, or cmd)."""
    workdir = tmpdir.mkdir("test_workdir")
    config_file = workdir.join("config.yaml")
    
    # Create minimal config
    test_config = {"model_args": {"zdim": 4}}
    utils.save_yaml(test_config, str(config_file))
    
    args = argparse.Namespace(workdir=str(workdir))
    
    # Should not raise any errors
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        view_config.main(args)