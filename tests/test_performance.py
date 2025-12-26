"""Performance tests and benchmarks for cryoDRGN.

This module provides performance testing and benchmarking for key cryoDRGN operations
to detect performance regressions and establish performance baselines.
"""

import pytest
import time
import numpy as np
import torch
from pathlib import Path
import psutil
import os
from typing import Dict, Any, List

from cryodrgn import dataset, ctf, fft, utils
from cryodrgn.source import ImageSource
from cryodrgn.lattice import Lattice
from tests.test_config import (
    isolated_temp_dir, create_mock_particles_file, 
    create_mock_poses_file, create_mock_ctf_file,
    PERFORMANCE_ITERATIONS, PERFORMANCE_WARMUP
)


class PerformanceBenchmark:
    """Helper class for performance benchmarking."""
    
    def __init__(self, name: str):
        self.name = name
        self.times = []
        self.memory_usage = []
        self.start_time = None
        self.process = psutil.Process()
    
    def __enter__(self):
        """Start timing and memory monitoring."""
        # Warm up
        torch.cuda.synchronize() if torch.cuda.is_available() else None
        
        self.start_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        self.start_time = time.perf_counter()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop timing and memory monitoring."""
        torch.cuda.synchronize() if torch.cuda.is_available() else None
        
        end_time = time.perf_counter()
        end_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        
        execution_time = end_time - self.start_time
        memory_used = end_memory - self.start_memory
        
        self.times.append(execution_time)
        self.memory_usage.append(memory_used)
    
    def summary(self) -> Dict[str, float]:
        """Get performance summary statistics."""
        return {
            "name": self.name,
            "mean_time": np.mean(self.times),
            "std_time": np.std(self.times),
            "min_time": np.min(self.times),
            "max_time": np.max(self.times),
            "mean_memory": np.mean(self.memory_usage),
            "max_memory": np.max(self.memory_usage),
        }


@pytest.mark.performance
@pytest.mark.benchmark
class TestDatasetPerformance:
    """Performance tests for dataset operations."""
    
    @pytest.mark.parametrize("n_particles,box_size", [
        (100, 64),
        (500, 128),
        (1000, 256),
    ])
    def test_dataset_loading_performance(self, isolated_temp_dir, 
                                       n_particles: int, box_size: int, 
                                       benchmark):
        """Benchmark dataset loading performance."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs",
            n_particles=n_particles,
            box_size=box_size
        )
        
        def load_dataset():
            return dataset.ImageDataset(str(particles_file), lazy=False)
        
        # Use pytest-benchmark if available, otherwise manual timing
        if benchmark:
            result = benchmark(load_dataset)
            assert result.N == n_particles
        else:
            bench = PerformanceBenchmark(f"dataset_loading_{n_particles}_{box_size}")
            for _ in range(PERFORMANCE_ITERATIONS):
                with bench:
                    data = load_dataset()
                    assert data.N == n_particles
            
            summary = bench.summary()
            # Basic performance assertions - adjust thresholds as needed
            assert summary["mean_time"] < 10.0, f"Dataset loading too slow: {summary}"
    
    @pytest.mark.parametrize("batch_size", [1, 8, 32, 128])
    def test_dataset_batch_access_performance(self, isolated_temp_dir, 
                                            batch_size: int, benchmark):
        """Benchmark dataset batch access performance."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs",
            n_particles=1000,
            box_size=64
        )
        
        data = dataset.ImageDataset(str(particles_file))
        indices = np.random.randint(0, 1000, batch_size)
        
        def access_batch():
            particles, _, _ = data[indices]
            return particles.shape[0]
        
        if benchmark:
            result = benchmark(access_batch)
            assert result == batch_size
        else:
            bench = PerformanceBenchmark(f"batch_access_{batch_size}")
            for _ in range(PERFORMANCE_ITERATIONS):
                with bench:
                    count = access_batch()
                    assert count == batch_size
            
            summary = bench.summary()
            # Performance should scale reasonably with batch size
            time_per_particle = summary["mean_time"] / batch_size
            assert time_per_particle < 0.01, f"Batch access too slow: {summary}"


@pytest.mark.performance
class TestFFTPerformance:
    """Performance tests for FFT operations."""
    
    @pytest.mark.parametrize("size", [64, 128, 256, 512])
    def test_fft_performance(self, size: int, benchmark):
        """Benchmark FFT operations."""
        # Create test data
        data = torch.randn(100, size, size, dtype=torch.float32)
        if torch.cuda.is_available():
            data = data.cuda()
        
        def run_fft():
            result = fft.fft2_center(data)
            return result.shape
        
        if benchmark:
            result = benchmark(run_fft)
            assert result == (100, size, size)
        else:
            bench = PerformanceBenchmark(f"fft_{size}")
            for _ in range(PERFORMANCE_ITERATIONS):
                with bench:
                    shape = run_fft()
                    assert shape == (100, size, size)
            
            summary = bench.summary()
            # FFT should be reasonably fast
            time_per_image = summary["mean_time"] / 100
            assert time_per_image < 0.1, f"FFT too slow: {summary}"
    
    @pytest.mark.parametrize("size", [64, 128, 256])
    def test_ifft_performance(self, size: int, benchmark):
        """Benchmark inverse FFT operations."""
        # Create test data in frequency domain
        data = torch.randn(50, size, size, dtype=torch.complex64)
        if torch.cuda.is_available():
            data = data.cuda()
        
        def run_ifft():
            result = fft.ifft2_center(data)
            return result.shape
        
        if benchmark:
            result = benchmark(run_ifft)
            assert result == (50, size, size)
        else:
            bench = PerformanceBenchmark(f"ifft_{size}")
            for _ in range(PERFORMANCE_ITERATIONS):
                with bench:
                    shape = run_ifft()
                    assert shape == (50, size, size)
            
            summary = bench.summary()
            time_per_image = summary["mean_time"] / 50
            assert time_per_image < 0.1, f"IFFT too slow: {summary}"


@pytest.mark.performance
class TestCTFPerformance:
    """Performance tests for CTF operations."""
    
    @pytest.mark.parametrize("n_images,box_size", [
        (100, 64),
        (500, 128),
        (100, 256),
    ])
    def test_ctf_computation_performance(self, isolated_temp_dir,
                                       n_images: int, box_size: int, 
                                       benchmark):
        """Benchmark CTF computation performance."""
        # Create lattice and CTF parameters
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        lattice = Lattice(box_size, device=device)
        
        # Create CTF parameters
        ctf_params = torch.tensor([
            [1.0,  # Apix
             np.random.uniform(1.0, 4.0, n_images),  # defocus_u
             np.random.uniform(1.0, 4.0, n_images),  # defocus_v
             np.random.uniform(0, 180, n_images),    # defocus_angle
             np.full(n_images, 300.0),               # voltage
             np.full(n_images, 2.7),                 # cs
             np.full(n_images, 0.1),                 # w
             np.zeros(n_images)                      # phase_shift
            ]
        ]).T.to(device)
        
        freqs = lattice.freqs2d.unsqueeze(0).expand(n_images, -1, -1)
        
        def compute_ctf_batch():
            results = []
            for i in range(n_images):
                c = ctf.compute_ctf(
                    freqs[i] / ctf_params[i, 0],
                    ctf_params[i, 1],
                    ctf_params[i, 2],
                    ctf_params[i, 3],
                    ctf_params[i, 4],
                    ctf_params[i, 5],
                    ctf_params[i, 6]
                )
                results.append(c)
            return len(results)
        
        if benchmark:
            result = benchmark(compute_ctf_batch)
            assert result == n_images
        else:
            bench = PerformanceBenchmark(f"ctf_computation_{n_images}_{box_size}")
            for _ in range(PERFORMANCE_ITERATIONS):
                with bench:
                    count = compute_ctf_batch()
                    assert count == n_images
            
            summary = bench.summary()
            time_per_image = summary["mean_time"] / n_images
            assert time_per_image < 0.01, f"CTF computation too slow: {summary}"


@pytest.mark.performance
class TestImageSourcePerformance:
    """Performance tests for ImageSource operations."""
    
    @pytest.mark.parametrize("n_particles", [100, 1000, 5000])
    def test_image_source_loading_performance(self, isolated_temp_dir, 
                                            n_particles: int, benchmark):
        """Benchmark ImageSource loading performance."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs",
            n_particles=n_particles,
            box_size=64
        )
        
        def load_image_source():
            src = ImageSource.from_file(str(particles_file), lazy=True)
            return src.n
        
        if benchmark:
            result = benchmark(load_image_source)
            assert result == n_particles
        else:
            bench = PerformanceBenchmark(f"image_source_loading_{n_particles}")
            for _ in range(PERFORMANCE_ITERATIONS):
                with bench:
                    count = load_image_source()
                    assert count == n_particles
            
            summary = bench.summary()
            # Loading should be fast regardless of dataset size (lazy loading)
            assert summary["mean_time"] < 5.0, f"ImageSource loading too slow: {summary}"
    
    @pytest.mark.parametrize("chunk_size", [10, 50, 100, 500])
    def test_chunked_access_performance(self, isolated_temp_dir, 
                                      chunk_size: int, benchmark):
        """Benchmark chunked access performance."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs",
            n_particles=1000,
            box_size=64
        )
        
        src = ImageSource.from_file(str(particles_file), lazy=True)
        
        def access_chunks():
            total_images = 0
            for indices, chunk in src.chunks(chunksize=chunk_size):
                total_images += len(indices)
                if total_images >= 200:  # Process first 200 images
                    break
            return total_images
        
        if benchmark:
            result = benchmark(access_chunks)
            assert result >= 200
        else:
            bench = PerformanceBenchmark(f"chunked_access_{chunk_size}")
            for _ in range(PERFORMANCE_ITERATIONS):
                with bench:
                    count = access_chunks()
                    assert count >= 200
            
            summary = bench.summary()
            time_per_image = summary["mean_time"] / 200
            assert time_per_image < 0.05, f"Chunked access too slow: {summary}"


@pytest.mark.performance 
@pytest.mark.slow
class TestMemoryUsage:
    """Test memory usage patterns and detect memory leaks."""
    
    def test_dataset_memory_usage(self, isolated_temp_dir):
        """Test memory usage of dataset operations."""
        particles_file = create_mock_particles_file(
            isolated_temp_dir / "particles.mrcs",
            n_particles=500,
            box_size=128
        )
        
        process = psutil.Process()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Load dataset multiple times and check for memory leaks
        for i in range(10):
            data = dataset.ImageDataset(str(particles_file), lazy=False)
            # Access some data
            particles, _, _ = data[np.random.randint(0, 500, 10)]
            del data, particles  # Explicit cleanup
            
            current_memory = process.memory_info().rss / 1024 / 1024  # MB
            memory_growth = current_memory - initial_memory
            
            # Memory should not grow indefinitely
            assert memory_growth < 1000, f"Memory leak detected: {memory_growth:.1f} MB growth"
    
    def test_fft_memory_usage(self):
        """Test memory usage of FFT operations."""
        if not torch.cuda.is_available():
            pytest.skip("GPU not available for memory test")
        
        device = torch.device('cuda')
        initial_memory = torch.cuda.memory_allocated(device)
        
        # Perform multiple FFT operations
        for i in range(10):
            data = torch.randn(100, 256, 256, device=device)
            fft_data = fft.fft2_center(data)
            ifft_data = fft.ifft2_center(fft_data)
            del data, fft_data, ifft_data
            torch.cuda.empty_cache()
        
        final_memory = torch.cuda.memory_allocated(device)
        memory_growth = (final_memory - initial_memory) / 1024 / 1024  # MB
        
        # Should not have significant memory growth
        assert memory_growth < 100, f"GPU memory leak: {memory_growth:.1f} MB growth"


def run_performance_suite(output_dir: Path = None):
    """Run the full performance test suite and generate a report."""
    if output_dir is None:
        output_dir = Path("performance_results")
    output_dir.mkdir(exist_ok=True)
    
    # Run performance tests and collect results
    import subprocess
    import json
    
    cmd = [
        "python", "-m", "pytest", 
        "tests/test_performance.py::TestDatasetPerformance",
        "tests/test_performance.py::TestFFTPerformance", 
        "tests/test_performance.py::TestCTFPerformance",
        "--benchmark-json", str(output_dir / "benchmark_results.json"),
        "-v"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Save detailed results
    with open(output_dir / "performance_log.txt", "w") as f:
        f.write("Performance Test Results\n")
        f.write("=" * 50 + "\n")
        f.write(f"Exit code: {result.returncode}\n\n")
        f.write("STDOUT:\n" + result.stdout + "\n\n")
        f.write("STDERR:\n" + result.stderr + "\n")
    
    print(f"Performance results saved to {output_dir}")
    return result.returncode == 0


if __name__ == "__main__":
    """Run performance tests when called as a script."""
    success = run_performance_suite()
    exit(0 if success else 1)