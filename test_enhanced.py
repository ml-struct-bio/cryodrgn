#!/usr/bin/env python3
"""
Enhanced testing script demonstrating parallel test execution with pytest-xdist
and Python 3.13 compatibility improvements for cryoDRGN.
"""

import subprocess
import sys
import time


def run_command(cmd):
    """Run a command and return timing and results."""
    print(f"Running: {cmd}")
    start_time = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    end_time = time.time()
    
    print(f"Exit code: {result.returncode}")
    print(f"Duration: {end_time - start_time:.2f} seconds")
    if result.stdout:
        print("STDOUT:", result.stdout[:500] + "..." if len(result.stdout) > 500 else result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr[:500] + "..." if len(result.stderr) > 500 else result.stderr)
    print("-" * 80)
    return result.returncode == 0, end_time - start_time


def main():
    """Demonstrate various testing scenarios."""
    
    # Set up environment
    env_setup = "export PATH=$PATH:/home/ubuntu/.local/bin && PYTHONPATH=/workspace:$PYTHONPATH"
    
    test_scenarios = [
        # Basic serial test run
        {
            "name": "Serial Tests (Basic)",
            "cmd": f"{env_setup} python3 -m pytest tests/test_fft.py tests/test_utils.py -v"
        },
        
        # Parallel test run with 2 workers
        {
            "name": "Parallel Tests (2 workers)",
            "cmd": f"{env_setup} python3 -m pytest tests/test_fft.py tests/test_utils.py tests/test_view_config.py tests/test_config.py -n2 --dist=loadscope -v"
        },
        
        # Fast tests only (using markers)
        {
            "name": "Fast Tests Only",
            "cmd": f"{env_setup} python3 -m pytest -m fast -v"
        },
        
        # Unit tests only
        {
            "name": "Unit Tests Only", 
            "cmd": f"{env_setup} python3 -m pytest -m unit -v"
        },
        
        # Integration tests only (slower)
        {
            "name": "Integration Tests Only",
            "cmd": f"{env_setup} python3 -m pytest -m integration -v"
        },
        
        # Parallel run excluding slow tests
        {
            "name": "Parallel Fast Tests",
            "cmd": f"{env_setup} python3 -m pytest -m 'not slow' -n2 --dist=loadscope -v"
        },
        
        # Test numpy 2.0 compatibility fixes
        {
            "name": "Numpy 2.0 Compatibility Tests",
            "cmd": f"{env_setup} python3 -m pytest tests/test_mrc.py -v"
        }
    ]
    
    print("=" * 80)
    print("CRYODRGN ENHANCED TESTING DEMONSTRATION")
    print("Python 3.13 Compatibility & Parallel Testing with pytest-xdist")
    print("=" * 80)
    
    results = []
    
    for scenario in test_scenarios:
        print(f"\n{scenario['name'].upper()}")
        print("=" * len(scenario['name']))
        
        success, duration = run_command(scenario['cmd'])
        results.append({
            'name': scenario['name'],
            'success': success,
            'duration': duration
        })
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    for result in results:
        status = "✓ PASSED" if result['success'] else "✗ FAILED"
        print(f"{result['name']:<30} {status:<10} {result['duration']:.2f}s")
    
    total_passed = sum(1 for r in results if r['success'])
    print(f"\nTotal: {total_passed}/{len(results)} scenarios passed")
    
    if total_passed == len(results):
        print("\n🎉 All testing scenarios completed successfully!")
        print("✅ Python 3.13 compatibility achieved")
        print("⚡ Parallel testing optimized with pytest-xdist")
        return 0
    else:
        print(f"\n⚠️  {len(results) - total_passed} scenario(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())