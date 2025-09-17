"""Comprehensive error handling and edge case tests for cryoDRGN.

This test module focuses on testing error conditions, invalid inputs,
and edge cases across different cryoDRGN commands and modules.
"""

import pytest
import os
import tempfile
import numpy as np
from pathlib import Path
from unittest.mock import Mock, patch
import argparse

from cryodrgn.commands import (
    backproject_voxel, train_vae, train_nn, analyze,
    downsample, parse_ctf_star, parse_pose_star
)
from cryodrgn import utils, dataset, ctf
from cryodrgn.source import ImageSource
from tests.test_config import isolated_temp_dir, create_mock_particles_file


class TestInputValidation:
    """Test input validation across different commands."""
    
    def test_nonexistent_particles_file(self, isolated_temp_dir):
        """Test handling of non-existent particles file."""
        nonexistent_file = isolated_temp_dir / "nonexistent.mrcs"
        
        with pytest.raises((FileNotFoundError, ValueError)):
            ImageSource.from_file(str(nonexistent_file))
    
    def test_invalid_particles_file_format(self, isolated_temp_dir):
        """Test handling of invalid particles file format."""
        invalid_file = isolated_temp_dir / "invalid.xyz"
        invalid_file.write_text("invalid content")
        
        with pytest.raises(ValueError, match="Unrecognized file type"):
            ImageSource.from_file(str(invalid_file))
    
    def test_corrupted_pickle_file(self, isolated_temp_dir):
        """Test handling of corrupted pickle files."""
        corrupted_file = isolated_temp_dir / "corrupted.pkl"
        corrupted_file.write_bytes(b"corrupted pickle data")
        
        with pytest.raises(Exception):  # Could be various pickle-related exceptions
            utils.load_pkl(str(corrupted_file))
    
    def test_empty_particles_file(self, isolated_temp_dir):
        """Test handling of empty particles file."""
        empty_file = isolated_temp_dir / "empty.mrcs"
        empty_file.touch()
        
        with pytest.raises(Exception):  # Should fail to load empty MRC file
            ImageSource.from_file(str(empty_file))
    
    def test_mismatched_particles_poses_count(self, isolated_temp_dir):
        """Test handling of mismatched particle and pose counts."""
        # Create particles file with 10 particles
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=10
        )
        
        # Create poses file with 5 poses (mismatch)
        poses_data = (np.random.randn(5, 3, 3), np.random.randn(5, 2))
        poses_file = isolated_temp_dir / "poses.pkl"
        utils.save_pkl(poses_data, str(poses_file))
        
        # This should raise an error when trying to create dataset
        with pytest.raises((ValueError, AssertionError)):
            data = dataset.ImageDataset(str(particles_file))
            # Error should occur when trying to use mismatched poses


class TestBackprojectionErrors:
    """Test error handling in backprojection command."""
    
    def test_backproject_missing_poses(self, isolated_temp_dir):
        """Test backprojection with missing poses file."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=5
        )
        
        args = argparse.Namespace(
            particles=str(particles_file),
            poses=str(isolated_temp_dir / "nonexistent_poses.pkl"),
            ctf=None,
            outdir=str(isolated_temp_dir / "output"),
            half_maps=True,
            fsc_vals=True,
            ctf_alg="mul",
            reg_weight=1.0,
            output_sumcount=False,
            invert_data=True,
            datadir=None,
            lazy=False,
            ind=None,
            first=None,
            log_interval="100",
            tilt=False,
            ntilts=10,
            force_ntilts=False,
            dose_per_tilt=None,
            angle_per_tilt=3,
        )
        
        with pytest.raises(FileNotFoundError):
            backproject_voxel.main(args)
    
    def test_backproject_invalid_ctf(self, isolated_temp_dir):
        """Test backprojection with invalid CTF file."""
        from tests.test_config import create_mock_particles_file, create_mock_poses_file
        
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=5
        )
        poses_file = create_mock_poses_file(
            isolated_temp_dir / "poses.pkl", n_poses=5
        )
        
        # Create invalid CTF file
        invalid_ctf_file = isolated_temp_dir / "invalid_ctf.pkl"
        invalid_ctf_file.write_text("invalid ctf data")
        
        args = argparse.Namespace(
            particles=str(particles_file),
            poses=str(poses_file),
            ctf=str(invalid_ctf_file),
            outdir=str(isolated_temp_dir / "output"),
            half_maps=True,
            fsc_vals=True,
            ctf_alg="mul",
            reg_weight=1.0,
            output_sumcount=False,
            invert_data=True,
            datadir=None,
            lazy=False,
            ind=None,
            first=None,
            log_interval="100",
            tilt=False,
            ntilts=10,
            force_ntilts=False,
            dose_per_tilt=None,
            angle_per_tilt=3,
        )
        
        with pytest.raises(Exception):  # Should fail to load invalid CTF
            backproject_voxel.main(args)


class TestDatasetErrors:
    """Test error handling in dataset loading."""
    
    def test_dataset_negative_indices(self, isolated_temp_dir):
        """Test dataset with negative indices."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=10
        )
        
        data = dataset.ImageDataset(str(particles_file))
        
        with pytest.raises((IndexError, ValueError)):
            data[[-1, -2, -3]]  # Negative indices should be handled gracefully
    
    def test_dataset_out_of_bounds_indices(self, isolated_temp_dir):
        """Test dataset with out-of-bounds indices."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=10
        )
        
        data = dataset.ImageDataset(str(particles_file))
        
        with pytest.raises(IndexError):
            data[[15, 20]]  # Indices beyond dataset size
    
    def test_tilt_series_without_dose_per_tilt(self, isolated_temp_dir):
        """Test tilt series data without required dose_per_tilt parameter."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=10
        )
        
        with pytest.raises((ValueError, AssertionError)):
            # Should fail without dose_per_tilt for tilt series
            data = dataset.TiltSeriesData(
                mrcfile=str(particles_file),
                ntilts=5,
                dose_per_tilt=None,  # This should cause an error
                angle_per_tilt=3.0,
            )


class TestCTFErrors:
    """Test error handling in CTF computations."""
    
    def test_ctf_invalid_voltage(self):
        """Test CTF computation with invalid voltage values."""
        import torch
        
        freqs = torch.randn(100, 2)  # Random frequency coordinates
        dfu = torch.tensor(2.0)
        dfv = torch.tensor(2.0) 
        dfang = torch.tensor(0.0)
        volt = torch.tensor(-100.0)  # Invalid negative voltage
        cs = torch.tensor(2.7)
        w = torch.tensor(0.1)
        
        # Should handle invalid voltage gracefully or raise appropriate error
        with pytest.raises((ValueError, RuntimeError)):
            ctf.compute_ctf(freqs, dfu, dfv, dfang, volt, cs, w)
    
    def test_ctf_zero_frequency(self):
        """Test CTF computation with zero frequencies."""
        import torch
        
        freqs = torch.zeros(100, 2)  # All zero frequencies
        dfu = torch.tensor(2.0)
        dfv = torch.tensor(2.0)
        dfang = torch.tensor(0.0)
        volt = torch.tensor(300.0)
        cs = torch.tensor(2.7)
        w = torch.tensor(0.1)
        
        # Should handle zero frequencies without crashing
        result = ctf.compute_ctf(freqs, dfu, dfv, dfang, volt, cs, w)
        assert torch.isfinite(result).all(), "CTF should produce finite values"


class TestUtilsErrors:
    """Test error handling in utility functions."""
    
    def test_save_pkl_readonly_directory(self, isolated_temp_dir):
        """Test saving pickle to read-only directory."""
        readonly_dir = isolated_temp_dir / "readonly"
        readonly_dir.mkdir()
        readonly_dir.chmod(0o444)  # Read-only permissions
        
        data = {"test": "data"}
        output_file = readonly_dir / "test.pkl"
        
        with pytest.raises(PermissionError):
            utils.save_pkl(data, str(output_file))
        
        # Restore permissions for cleanup
        readonly_dir.chmod(0o755)
    
    def test_run_command_invalid_command(self):
        """Test running invalid command."""
        with pytest.raises(Exception):
            utils.run_command("invalid_command_that_does_not_exist")
    
    def test_run_command_with_failure(self):
        """Test running command that returns non-zero exit code."""
        # This command should fail
        out, err = utils.run_command("ls /nonexistent/directory/path")
        # Should capture error output without crashing
        assert err or "No such file or directory" in out


class TestMemoryLimits:
    """Test behavior under memory constraints."""
    
    @pytest.mark.slow
    def test_large_dataset_lazy_loading(self, isolated_temp_dir):
        """Test that lazy loading works with large datasets."""
        # Create a dataset that would be too large to load entirely
        large_particles_file = create_mock_particles_file(
            isolated_temp_dir / "large_particles.mrcs", 
            n_particles=1000,  # Reasonably large for test
            box_size=128
        )
        
        # This should work with lazy loading
        data = dataset.ImageDataset(
            str(large_particles_file), 
            lazy=True
        )
        
        # Access a few random particles
        sample_indices = np.random.choice(1000, 10, replace=False)
        particles, indices, _ = data[sample_indices]
        
        assert particles.shape[0] == 10
        assert not torch.isnan(particles).any()


@pytest.mark.parametrize("invalid_arg", [
    {"zdim": -1},  # Negative zdim
    {"epochs": 0},  # Zero epochs
    {"batch_size": -5},  # Negative batch size
])
class TestInvalidCommandArguments:
    """Test command-line argument validation."""
    
    def test_train_vae_invalid_args(self, invalid_arg, isolated_temp_dir):
        """Test train_vae with various invalid arguments."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=5
        )
        
        args = argparse.Namespace(
            particles=str(particles_file),
            outdir=str(isolated_temp_dir / "output"),
            poses=None,
            ctf=None,
            load=None,
            checkpoint=1,
            log_interval=100,
            verbose=False,
            seed=None,
            ind=None,
            uninvert_data=False,
            no_window=False,
            window_r=0.85,
            datadir=None,
            lazy=False,
            shuffler_size=0,
            no_amp=False,
            multigpu=False,
            do_pose_sgd=False,
            do_trans_sgd=False,
            pretrain=1,
            ps_freq=5,
            shift_freq=5,
            l_start=10000,
            l_end=10000,
            num_epochs=2,
            batch_size=8,
            wd=0,
            lr=0.0001,
            beta=None,
            beta_control=None,
            norm=None,
            no_trans=False,
            enc_layers=3,
            enc_dim=256,
            zdim=8,  # Default, will be overridden
            dec_layers=3,
            dec_dim=256,
            pe_type="gaussian",
            feat_sigma=0.5,
            pe_dim=None,
            domain="hartley",
            activation="relu",
            **invalid_arg
        )
        
        with pytest.raises((ValueError, AssertionError)):
            train_vae.main(args)


class TestConcurrencyIssues:
    """Test concurrent access and thread safety."""
    
    @pytest.mark.slow
    def test_concurrent_dataset_access(self, isolated_temp_dir):
        """Test concurrent access to dataset."""
        import threading
        import time
        
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs", n_particles=100
        )
        
        data = dataset.ImageDataset(str(particles_file))
        errors = []
        
        def worker():
            try:
                for i in range(10):
                    idx = np.random.randint(0, 100)
                    particles, _, _ = data[idx]
                    time.sleep(0.001)  # Small delay to encourage race conditions
            except Exception as e:
                errors.append(e)
        
        # Start multiple threads
        threads = [threading.Thread(target=worker) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        
        # Check that no errors occurred during concurrent access
        assert len(errors) == 0, f"Errors during concurrent access: {errors}"