#!/usr/bin/env python3
"""
Optimized test runner for cryodrgn_beta with improved parallelization strategies.

This script provides different test execution strategies to maximize throughput
on systems with 4 CPUs while managing resource usage effectively.
"""

import subprocess
import sys
import time
import argparse
from pathlib import Path


def run_command(cmd, cwd=None):
    """Run a command and return execution time and result."""
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=cwd, capture_output=True, text=True, check=True
        )
        duration = time.time() - start_time
        return True, duration, result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        duration = time.time() - start_time
        return False, duration, e.stdout, e.stderr


def run_stratified_tests():
    """
    Run tests in a stratified manner - separate fast and slow tests
    for optimal parallelization.
    """
    print("🚀 Running Stratified Test Strategy")
    print("=" * 60)
    
    # Strategy 1: Run fast unit tests first with high parallelization
    print("\n📝 Phase 1: Fast Unit Tests (high parallelization)")
    fast_tests = [
        "test_utils.py",
        "test_fft.py", 
        "test_source.py",
        "test_mrc.py",
        "test_masks.py",
        "test_flip_hand.py",
        "test_invert_contrast.py",
        "test_phase_flip.py",
        "test_entropy.py",
        "test_view_*.py",
        "test_fsc.py",
        "test_pc_traversal.py",
        "test_direct_traversal.py"
    ]
    
    fast_cmd = f"python3 -m pytest -n4 --dist=each -v " + " ".join([f"tests/{test}" for test in fast_tests])
    success, duration, stdout, stderr = run_command(fast_cmd)
    print(f"⏱️  Fast tests completed in {duration:.2f}s")
    if not success:
        print(f"❌ Fast tests failed:\n{stderr}")
        return False
    
    # Strategy 2: Run medium complexity tests with moderate parallelization  
    print("\n⚙️  Phase 2: Medium Complexity Tests (moderate parallelization)")
    medium_tests = [
        "test_parse.py",
        "test_relion.py", 
        "test_writestar.py",
        "test_dataset.py",
        "test_downsample.py",
        "test_translate.py",
        "test_clean.py",
        "test_filter_*.py",
        "test_select_*.py",
        "test_add_psize.py",
        "test_graph_traversal.py",
        "test_eval_images.py"
    ]
    
    medium_cmd = f"python3 -m pytest -n3 --dist=loadscope -v " + " ".join([f"tests/{test}" for test in medium_tests])
    success, duration, stdout, stderr = run_command(medium_cmd)
    print(f"⏱️  Medium tests completed in {duration:.2f}s")
    if not success:
        print(f"❌ Medium tests failed:\n{stderr}")
        return False
    
    # Strategy 3: Run heavy computational tests with lower parallelization
    print("\n🏋️  Phase 3: Heavy Computational Tests (conservative parallelization)")
    heavy_tests = [
        "test_integration.py",
        "test_reconstruct_*.py", 
        "test_backprojection.py",
        "test_read_filter_write.py"
    ]
    
    heavy_cmd = f"python3 -m pytest -n2 --dist=loadgroup -v " + " ".join([f"tests/{test}" for test in heavy_tests])
    success, duration, stdout, stderr = run_command(heavy_cmd)
    print(f"⏱️  Heavy tests completed in {duration:.2f}s")
    if not success:
        print(f"❌ Heavy tests failed:\n{stderr}")
        return False
    
    print("\n✅ All test phases completed successfully!")
    return True


def run_optimized_single_pass():
    """
    Run all tests in a single pass with optimized settings for 4 CPUs.
    """
    print("🚀 Running Optimized Single Pass Strategy")
    print("=" * 60)
    
    # Use 4 workers with 'each' distribution for better load balancing
    cmd = "python3 -m pytest -n4 --dist=each -v --tb=short --maxfail=10"
    
    success, duration, stdout, stderr = run_command(cmd)
    print(f"⏱️  All tests completed in {duration:.2f}s")
    
    if success:
        print("✅ All tests passed!")
        print(stdout)
    else:
        print(f"❌ Some tests failed:\n{stderr}")
        
    return success


def run_coverage_analysis():
    """
    Run tests with coverage analysis to identify gaps.
    """
    print("📊 Running Coverage Analysis")
    print("=" * 60)
    
    # Install coverage if not available
    subprocess.run(["python3", "-m", "pip", "install", "--user", "coverage", "pytest-cov", "--break-system-packages"], 
                   capture_output=True)
    
    cmd = "python3 -m pytest --cov=cryodrgn --cov-report=html --cov-report=term-missing -n2 --dist=loadscope"
    
    success, duration, stdout, stderr = run_command(cmd)
    print(f"⏱️  Coverage analysis completed in {duration:.2f}s")
    
    if success:
        print("✅ Coverage analysis complete! Check htmlcov/ directory for detailed report.")
        print(stdout)
    else:
        print(f"❌ Coverage analysis failed:\n{stderr}")
        
    return success


def main():
    parser = argparse.ArgumentParser(description="Optimized test runner for cryodrgn_beta")
    parser.add_argument(
        "--strategy", 
        choices=["stratified", "single-pass", "coverage"],
        default="stratified",
        help="Test execution strategy (default: stratified)"
    )
    
    args = parser.parse_args()
    
    # Ensure we're in the right directory
    if not Path("tests").exists():
        print("❌ Tests directory not found. Please run from project root.")
        sys.exit(1)
    
    start_time = time.time()
    
    if args.strategy == "stratified":
        success = run_stratified_tests()
    elif args.strategy == "single-pass":
        success = run_optimized_single_pass()
    elif args.strategy == "coverage":
        success = run_coverage_analysis()
    
    total_time = time.time() - start_time
    print(f"\n🏁 Total execution time: {total_time:.2f}s")
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()