"""
Optimizations for conftest.py to improve test performance and parallelization.

This file contains improvements that can be applied to the existing conftest.py
to better support parallel test execution with 4 CPUs.
"""

import pytest
import os
from pathlib import Path
from typing import Dict, Any
import pickle
import tempfile
import shutil
from functools import lru_cache


# Optimized fixture management
class TestDataCache:
    """Singleton cache for expensive test data to avoid reloading."""
    _instance = None
    _cache = {}
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    @lru_cache(maxsize=32)
    def get_test_data(self, file_path: str):
        """Cache commonly used test data files."""
        if file_path not in self._cache:
            if file_path.endswith('.pkl'):
                with open(file_path, 'rb') as f:
                    self._cache[file_path] = pickle.load(f)
            # Add other file type handlers as needed
        return self._cache[file_path]


# Session-scoped fixtures for expensive setup
@pytest.fixture(scope="session")
def test_data_cache():
    """Provide cached test data across the session."""
    return TestDataCache()


@pytest.fixture(scope="session")
def temp_workspace(tmp_path_factory):
    """Create a shared temporary workspace for the session."""
    workspace = tmp_path_factory.mktemp("cryodrgn_test_workspace")
    yield workspace
    # Cleanup handled automatically by tmp_path_factory


@pytest.fixture(scope="session")
def common_test_files(temp_workspace):
    """Copy commonly used test files to workspace once per session."""
    data_dir = Path(__file__).parent / "tests" / "data"
    common_files = [
        "toy_projections.mrcs",
        "toy_angles.pkl", 
        "toy_rot_trans.pkl",
        "test_ctf.pkl",
        "hand.mrcs",
        "hand_rot.pkl"
    ]
    
    file_paths = {}
    for filename in common_files:
        src = data_dir / filename
        if src.exists():
            dst = temp_workspace / filename
            shutil.copy2(src, dst)
            file_paths[filename] = str(dst)
    
    return file_paths


# Optimized data fixtures with better scoping
@pytest.fixture(scope="class")  # Changed from function to class scope
def particles_cached(request, common_test_files, test_data_cache):
    """Optimized particles fixture with caching."""
    # Implementation similar to original but with caching
    pass


@pytest.fixture(scope="class")  # Changed from function to class scope  
def poses_cached(request, common_test_files, test_data_cache):
    """Optimized poses fixture with caching."""
    # Implementation similar to original but with caching
    pass


# Parallelization helpers
def pytest_configure(config):
    """Configure pytest for optimal parallel execution."""
    # Add custom markers
    config.addinivalue_line("markers", "slow: mark test as slow running")
    config.addinivalue_line("markers", "training: mark test as requiring training")
    config.addinivalue_line("markers", "integration: mark test as integration test")
    config.addinivalue_line("markers", "parallel_safe: mark test as safe for parallel execution")
    
    # Set environment variables for better parallel performance
    os.environ["NUMEXPR_MAX_THREADS"] = "1"  # Prevent oversubscription
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"


def pytest_collection_modifyitems(config, items):
    """Modify test collection for better parallel distribution."""
    # Mark slow tests
    slow_markers = ["integration", "training", "backproject"]
    
    for item in items:
        # Auto-mark slow tests
        if any(marker in item.nodeid.lower() for marker in slow_markers):
            item.add_marker(pytest.mark.slow)
        
        # Auto-mark parallel safe tests
        if not any(marker in item.nodeid.lower() for marker in ["train", "integration"]):
            item.add_marker(pytest.mark.parallel_safe)


@pytest.fixture(scope="session", autouse=True)
def setup_test_environment():
    """Set up optimal test environment once per session."""
    # Configure matplotlib for headless testing
    import matplotlib
    matplotlib.use('Agg')
    
    # Set random seeds for reproducibility
    import numpy as np
    import torch
    np.random.seed(42)
    torch.manual_seed(42)
    
    # Disable CUDA for consistent testing unless explicitly needed
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    
    yield
    
    # Cleanup
    pass


# Resource management for heavy tests
class ResourceManager:
    """Manage computational resources for heavy tests."""
    
    def __init__(self):
        self._heavy_test_lock = None
    
    @pytest.fixture(scope="session")
    def heavy_test_manager(self):
        """Manage heavy computational tests to prevent resource conflicts."""
        try:
            # Try to import threading for locks
            import threading
            self._heavy_test_lock = threading.Lock()
        except ImportError:
            self._heavy_test_lock = None
        
        return self
    
    def acquire_heavy_resource(self):
        """Acquire lock for heavy computational tests."""
        if self._heavy_test_lock:
            self._heavy_test_lock.acquire()
    
    def release_heavy_resource(self):
        """Release lock for heavy computational tests."""
        if self._heavy_test_lock:
            self._heavy_test_lock.release()


# Optimized TrainDir class for better parallel testing
class OptimizedTrainDir:
    """Optimized version of TrainDir for better parallel testing."""
    
    _cache = {}  # Class-level cache for training results
    
    def __init__(self, dataset: str, train_cmd: str, epochs: int = 5, **kwargs):
        self.cache_key = f"{dataset}_{train_cmd}_{epochs}"
        
        # Check if we have cached results
        if self.cache_key in self._cache:
            self.outdir = self._cache[self.cache_key]
            return
        
        # Otherwise create new training directory
        # Implementation here...
        # Store in cache for reuse
        self._cache[self.cache_key] = self.outdir