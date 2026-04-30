"""Test configuration and utilities for cryoDRGN tests.

This module provides common configuration, fixtures, and utilities
that are used across multiple test modules.
"""

import os
import pytest
import tempfile
import shutil
from typing import Generator, Any
from pathlib import Path

# Test configuration constants
TEST_DATA_DIR = Path(__file__).parent / "data"
TEMP_DIR_PREFIX = "cryodrgn_test_"
DEFAULT_TIMEOUT = 60  # seconds for slow tests

# Performance test configuration
PERFORMANCE_ITERATIONS = 3
PERFORMANCE_WARMUP = 1

# Common test parameters
SMALL_DATASET_SIZE = 5
MEDIUM_DATASET_SIZE = 50
LARGE_DATASET_SIZE = 1000

def skip_if_no_gpu():
    """Skip test if no GPU is available."""
    import torch
    return pytest.mark.skipif(
        not torch.cuda.is_available(),
        reason="GPU not available"
    )

def skip_if_slow(reason="Slow test"):
    """Skip test if SKIP_SLOW_TESTS environment variable is set."""
    return pytest.mark.skipif(
        os.environ.get("SKIP_SLOW_TESTS", "").lower() in ("1", "true", "yes"),
        reason=reason
    )

@pytest.fixture(scope="session")
def temp_workspace() -> Generator[Path, None, None]:
    """Create a temporary workspace for tests."""
    with tempfile.TemporaryDirectory(prefix=TEMP_DIR_PREFIX) as temp_dir:
        workspace = Path(temp_dir)
        yield workspace

@pytest.fixture(scope="function") 
def isolated_temp_dir() -> Generator[Path, None, None]:
    """Create an isolated temporary directory for each test function."""
    with tempfile.TemporaryDirectory(prefix=TEMP_DIR_PREFIX) as temp_dir:
        yield Path(temp_dir)

class TestMetrics:
    """Helper class for collecting test metrics and performance data."""
    
    def __init__(self):
        self.execution_times = []
        self.memory_usage = []
        self.errors = []
    
    def record_time(self, execution_time: float):
        """Record execution time for performance analysis."""
        self.execution_times.append(execution_time)
    
    def record_memory(self, memory_mb: float):
        """Record memory usage for performance analysis."""
        self.memory_usage.append(memory_mb)
    
    def record_error(self, error: Exception):
        """Record error for analysis."""
        self.errors.append(error)
    
    @property
    def avg_execution_time(self) -> float:
        """Get average execution time."""
        return sum(self.execution_times) / len(self.execution_times) if self.execution_times else 0.0
    
    @property
    def max_memory_usage(self) -> float:
        """Get maximum memory usage."""
        return max(self.memory_usage) if self.memory_usage else 0.0

def assert_files_equal(file1: Path, file2: Path, tolerance: float = 1e-6):
    """Assert that two files contain approximately equal numerical data."""
    import numpy as np
    
    # Handle different file types
    if file1.suffix == '.npy' and file2.suffix == '.npy':
        data1 = np.load(file1)
        data2 = np.load(file2)
        np.testing.assert_allclose(data1, data2, rtol=tolerance)
    elif file1.suffix in ['.pkl', '.pickle'] and file2.suffix in ['.pkl', '.pickle']:
        import pickle
        with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
            data1 = pickle.load(f1)
            data2 = pickle.load(f2)
        if isinstance(data1, np.ndarray) and isinstance(data2, np.ndarray):
            np.testing.assert_allclose(data1, data2, rtol=tolerance)
        else:
            assert data1 == data2
    else:
        # For text files, compare content
        assert file1.read_text() == file2.read_text()

def create_mock_particles_file(output_path: Path, n_particles: int = 10, 
                              box_size: int = 64) -> Path:
    """Create a mock particles file for testing."""
    import numpy as np
    from cryodrgn.source import write_mrc
    
    # Create random particle data
    particles_data = np.random.randn(n_particles, box_size, box_size).astype(np.float32)
    
    # Write to MRC format
    write_mrc(str(output_path), particles_data, Apix=1.0)
    return output_path

def create_mock_poses_file(output_path: Path, n_poses: int = 10) -> Path:
    """Create a mock poses file for testing."""
    import numpy as np
    import pickle
    
    # Create random rotation matrices and translations  
    from scipy.spatial.transform import Rotation
    rotations = Rotation.random(n_poses).as_matrix()
    translations = np.random.randn(n_poses, 2) * 5  # Random translations
    
    poses_data = (rotations, translations)
    
    with open(output_path, 'wb') as f:
        pickle.dump(poses_data, f)
    
    return output_path

def create_mock_ctf_file(output_path: Path, n_particles: int = 10) -> Path:
    """Create a mock CTF parameters file for testing."""
    import numpy as np
    import pickle
    
    # Create realistic CTF parameters
    # Format: [Apix, defocus_u, defocus_v, defocus_angle, voltage, cs, w, phase_shift]
    ctf_params = np.array([
        [1.0,  # Apix
         np.random.uniform(1.0, 4.0, n_particles),  # defocus_u
         np.random.uniform(1.0, 4.0, n_particles),  # defocus_v  
         np.random.uniform(0, 180, n_particles),    # defocus_angle
         np.full(n_particles, 300.0),               # voltage
         np.full(n_particles, 2.7),                 # cs
         np.full(n_particles, 0.1),                 # w
         np.zeros(n_particles)                      # phase_shift
        ]
    ]).T
    
    with open(output_path, 'wb') as f:
        pickle.dump(ctf_params, f)
    
    return output_path