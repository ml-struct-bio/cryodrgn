"""Parallel-optimized configuration for pytest execution.

This file provides improved fixtures and configurations for better
parallel test execution with pytest-xdist.
"""

import pytest
import tempfile
import os
import shutil
from pathlib import Path
import multiprocessing as mp
from unittest.mock import patch


@pytest.fixture(scope="session")
def parallel_tmp_dir():
    """Create a session-scoped temporary directory for parallel execution.
    
    This ensures each test session gets its own isolated temp directory.
    """
    with tempfile.TemporaryDirectory(prefix="cryodrgn_test_") as tmpdir:
        yield tmpdir


@pytest.fixture(scope="function")
def isolated_tmp_dir():
    """Create function-scoped temporary directory for test isolation.
    
    Each test function gets its own clean temporary directory.
    """
    with tempfile.TemporaryDirectory(prefix="cryodrgn_func_") as tmpdir:
        # Change to the temp directory for the test
        old_cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            yield tmpdir
        finally:
            os.chdir(old_cwd)


@pytest.fixture(scope="session")
def cpu_count():
    """Get the number of available CPU cores for parallel processing."""
    return mp.cpu_count()


@pytest.fixture(autouse=True)
def setup_test_environment():
    """Setup clean test environment for each test.
    
    This fixture automatically runs for every test to ensure clean state.
    """
    # Set environment variables for reproducible testing
    os.environ["PYTHONHASHSEED"] = "42"
    
    # Suppress matplotlib GUI
    os.environ["MPLBACKEND"] = "Agg"
    
    # Set torch to CPU mode for consistent testing
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    
    yield
    
    # Cleanup any environment changes if needed


@pytest.fixture(scope="session")
def test_data_dir():
    """Get the test data directory path."""
    return Path(__file__).parent / "data"


@pytest.fixture
def mock_gpu_available():
    """Mock GPU availability for testing GPU-related code paths."""
    with patch('torch.cuda.is_available', return_value=True):
        with patch('torch.cuda.device_count', return_value=2):
            yield


@pytest.fixture
def mock_no_gpu():
    """Mock no GPU availability for testing CPU-only code paths.""" 
    with patch('torch.cuda.is_available', return_value=False):
        with patch('torch.cuda.device_count', return_value=0):
            yield


# Pytest collection hooks for better parallel execution
def pytest_collection_modifyitems(config, items):
    """Modify collected test items to optimize parallel execution.
    
    This function runs during test collection and can be used to:
    - Mark tests for better parallel distribution
    - Skip tests based on conditions
    - Reorder tests for optimal execution
    """
    
    # Mark slow tests
    for item in items:
        # Mark integration tests as slow
        if "integration" in item.nodeid:
            item.add_marker(pytest.mark.slow)
        
        # Mark tests that involve model training as slow
        if any(keyword in item.nodeid.lower() for keyword in [
            "train", "abinit", "reconstruct", "notebook"
        ]):
            item.add_marker(pytest.mark.slow)
        
        # Mark parsing tests for better grouping
        if any(keyword in item.nodeid.lower() for keyword in [
            "parse", "ctf", "pose", "star", "csparc"
        ]):
            item.add_marker(pytest.mark.parsing)
        
        # Mark utility tests
        if any(keyword in item.nodeid.lower() for keyword in [
            "utils", "mrc", "fft", "mask", "fsc"
        ]):
            item.add_marker(pytest.mark.utilities)


def pytest_configure(config):
    """Configure pytest for better parallel execution."""
    
    # Register custom markers
    config.addinivalue_line("markers", "slow: marks tests as slow")
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line("markers", "unit: marks tests as unit tests")
    config.addinivalue_line("markers", "parsing: marks tests as parsing-related")
    config.addinivalue_line("markers", "utilities: marks tests as utility functions")
    config.addinivalue_line("markers", "reconstruction: marks tests as reconstruction-related")
    config.addinivalue_line("markers", "analysis: marks tests as analysis-related")
    config.addinivalue_line("markers", "forked: marks tests that should run in separate processes")
    
    # Set up logging for better debugging in parallel execution
    import logging
    logging.basicConfig(
        level=logging.WARNING,  # Reduce log noise in tests
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def pytest_runtest_setup(item):
    """Setup for individual test runs.
    
    This runs before each test and can be used for test-specific setup.
    """
    # Ensure each test starts with a clean state
    import gc
    gc.collect()


def pytest_runtest_teardown(item):
    """Teardown after individual test runs.
    
    This runs after each test for cleanup.
    """
    # Clean up any remaining resources
    import gc
    gc.collect()


# Session-scoped fixtures for expensive setup
@pytest.fixture(scope="session")
def shared_test_resources():
    """Create shared test resources that can be reused across tests.
    
    This fixture sets up expensive resources once per test session.
    """
    resources = {
        'temp_dir': tempfile.mkdtemp(prefix="cryodrgn_shared_"),
        'test_data_available': True
    }
    
    yield resources
    
    # Cleanup
    if os.path.exists(resources['temp_dir']):
        shutil.rmtree(resources['temp_dir'])


# Parallel execution optimization markers
pytest_plugins = [
    "xdist"  # Enable pytest-xdist plugin for parallel execution
]