#!/usr/bin/env python3
"""Comprehensive test runner for cryoDRGN testing suite.

This script provides a unified interface to run different categories of tests,
generate reports, and manage test environments.
"""

import os
import sys
import argparse
import subprocess
import time
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import json
import tempfile

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestRunner:
    """Main test runner class."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.test_dir = project_root / "tests"
        self.results_dir = project_root / "test_results"
        self.results_dir.mkdir(exist_ok=True)
        
        # Test categories and their corresponding files/markers
        self.test_categories = {
            "unit": {
                "description": "Fast unit tests",
                "markers": ["unit", "not slow", "not integration", "not performance"],
                "files": [
                    "test_dataset.py",
                    "test_fft.py", 
                    "test_utils.py",
                    "test_ctf.py",  # If exists
                    "test_source.py",
                ],
            },
            "integration": {
                "description": "Integration tests",
                "markers": ["integration"],
                "files": [
                    "test_integration.py",
                    "test_backprojection.py", 
                    "test_reconstruct_*.py",
                ],
            },
            "performance": {
                "description": "Performance and benchmark tests", 
                "markers": ["performance", "benchmark"],
                "files": ["test_performance.py"],
            },
            "error_handling": {
                "description": "Error handling and edge case tests",
                "markers": [],
                "files": ["test_error_handling.py"],
            },
            "property_based": {
                "description": "Property-based tests with Hypothesis",
                "markers": [],
                "files": ["test_property_based.py"],
            },
            "mocking": {
                "description": "Tests with mocked dependencies",
                "markers": [],
                "files": ["test_mocking.py"],
            },
            "documentation": {
                "description": "Documentation and docstring tests",
                "markers": [],
                "files": ["test_documentation.py"],
            },
            "slow": {
                "description": "Slow tests (long-running)",
                "markers": ["slow"],
                "files": [],
            },
            "gpu": {
                "description": "GPU-specific tests",
                "markers": ["gpu"],
                "files": [],
            },
            "cmdline": {
                "description": "Command-line interface tests",
                "markers": ["cmdline"],
                "files": [],
            },
        }
    
    def get_pytest_command(self, category: str, extra_args: List[str] = None) -> List[str]:
        """Generate pytest command for a given test category."""
        cmd = ["python", "-m", "pytest"]
        
        # Add coverage by default for most categories
        if category not in ["performance", "documentation"]:
            cmd.extend(["--cov=cryodrgn", "--cov-report=xml", "--cov-append"])
        
        # Add markers
        test_config = self.test_categories.get(category, {})
        markers = test_config.get("markers", [])
        if markers:
            marker_expr = " and ".join(markers)
            cmd.extend(["-m", marker_expr])
        
        # Add specific test files
        test_files = test_config.get("files", [])
        for pattern in test_files:
            if "*" in pattern:
                # Glob pattern - add all matching files
                matching_files = list(self.test_dir.glob(pattern))
                cmd.extend([str(f) for f in matching_files])
            else:
                # Specific file
                test_file = self.test_dir / pattern
                if test_file.exists():
                    cmd.append(str(test_file))
        
        # Add extra arguments
        if extra_args:
            cmd.extend(extra_args)
        
        # Default to test directory if no specific files
        if not any(str(self.test_dir) in arg or arg.endswith('.py') for arg in cmd):
            cmd.append(str(self.test_dir))
        
        # Add verbosity and output options
        cmd.extend([
            "-v",
            "--tb=short",
            "--color=yes",
            f"--junitxml={self.results_dir / f'{category}_results.xml'}",
        ])
        
        return cmd
    
    def run_test_category(self, category: str, extra_args: List[str] = None) -> Tuple[int, float, str]:
        """Run a specific test category and return (exit_code, duration, output)."""
        print(f"\n{'='*60}")
        print(f"Running {category} tests:")
        print(f"  {self.test_categories.get(category, {}).get('description', '')}")
        print(f"{'='*60}")
        
        start_time = time.time()
        
        cmd = self.get_pytest_command(category, extra_args)
        print(f"Command: {' '.join(cmd)}")
        
        # Run the tests
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=self.project_root
        )
        
        duration = time.time() - start_time
        
        # Save output
        output_file = self.results_dir / f"{category}_output.txt"
        with open(output_file, 'w') as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Exit code: {result.returncode}\n")
            f.write(f"Duration: {duration:.2f}s\n")
            f.write(f"\n--- STDOUT ---\n{result.stdout}")
            f.write(f"\n--- STDERR ---\n{result.stderr}")
        
        # Print results
        print(f"Exit code: {result.returncode}")
        print(f"Duration: {duration:.2f}s")
        print(f"Output saved to: {output_file}")
        
        if result.stdout:
            print("\nSTDOUT:")
            print(result.stdout[-1000:])  # Last 1000 chars
        
        if result.stderr and result.returncode != 0:
            print("\nSTDERR:")
            print(result.stderr[-500:])  # Last 500 chars
        
        return result.returncode, duration, result.stdout
    
    def run_all_categories(self, exclude: List[str] = None) -> Dict[str, Tuple[int, float]]:
        """Run all test categories and return results."""
        exclude = exclude or []
        results = {}
        
        print(f"Running comprehensive test suite")
        print(f"Excluding categories: {exclude}")
        
        for category in self.test_categories:
            if category in exclude:
                print(f"Skipping {category} tests (excluded)")
                continue
            
            exit_code, duration, _ = self.run_test_category(category)
            results[category] = (exit_code, duration)
        
        return results
    
    def generate_summary_report(self, results: Dict[str, Tuple[int, float]]):
        """Generate a summary report of test results."""
        print(f"\n{'='*60}")
        print("TEST SUMMARY REPORT")
        print(f"{'='*60}")
        
        total_duration = sum(duration for _, duration in results.values())
        total_passed = sum(1 for exit_code, _ in results.values() if exit_code == 0)
        total_failed = len(results) - total_passed
        
        print(f"Total test categories: {len(results)}")
        print(f"Passed: {total_passed}")
        print(f"Failed: {total_failed}")
        print(f"Total duration: {total_duration:.2f}s")
        print()
        
        # Detailed results
        print("Category Results:")
        print("-" * 60)
        for category, (exit_code, duration) in results.items():
            status = "PASS" if exit_code == 0 else "FAIL"
            description = self.test_categories.get(category, {}).get('description', '')
            print(f"{category:15} | {status:4} | {duration:6.2f}s | {description}")
        
        # Save summary to file
        summary_file = self.results_dir / "test_summary.json"
        summary_data = {
            "timestamp": time.time(),
            "total_categories": len(results),
            "passed": total_passed,
            "failed": total_failed,
            "total_duration": total_duration,
            "results": {
                category: {
                    "exit_code": exit_code,
                    "duration": duration,
                    "status": "pass" if exit_code == 0 else "fail"
                }
                for category, (exit_code, duration) in results.items()
            }
        }
        
        with open(summary_file, 'w') as f:
            json.dump(summary_data, f, indent=2)
        
        print(f"\nDetailed summary saved to: {summary_file}")
        
        return total_failed == 0
    
    def install_dependencies(self):
        """Install test dependencies."""
        print("Installing test dependencies...")
        cmd = [sys.executable, "-m", "pip", "install", "-e", ".[dev]"]
        result = subprocess.run(cmd, cwd=self.project_root)
        return result.returncode == 0
    
    def check_environment(self):
        """Check test environment and dependencies."""
        print("Checking test environment...")
        
        issues = []
        
        # Check Python version
        if sys.version_info < (3, 10):
            issues.append(f"Python version {sys.version_info} < 3.10")
        
        # Check required packages
        required_packages = [
            "pytest", "pytest-cov", "pytest-benchmark", "pytest-mock", 
            "hypothesis", "torch", "numpy", "scipy"
        ]
        
        for package in required_packages:
            try:
                __import__(package.replace('-', '_'))
            except ImportError:
                issues.append(f"Missing package: {package}")
        
        # Check GPU availability
        try:
            import torch
            gpu_available = torch.cuda.is_available()
            print(f"GPU available: {gpu_available}")
            if gpu_available:
                print(f"GPU device count: {torch.cuda.device_count()}")
                print(f"GPU device name: {torch.cuda.get_device_name(0)}")
        except:
            issues.append("Cannot check GPU availability")
        
        if issues:
            print("Environment issues found:")
            for issue in issues:
                print(f"  - {issue}")
            return False
        else:
            print("Environment check passed!")
            return True
    
    def quick_test(self) -> bool:
        """Run a quick subset of tests for development."""
        print("Running quick development tests...")
        
        quick_categories = ["unit", "error_handling"]
        results = {}
        
        for category in quick_categories:
            exit_code, duration, _ = self.run_test_category(
                category, 
                extra_args=["--maxfail=3", "-x"]  # Stop on first few failures
            )
            results[category] = (exit_code, duration)
        
        return self.generate_summary_report(results)


def main():
    """Main entry point for test runner."""
    parser = argparse.ArgumentParser(
        description="cryoDRGN Comprehensive Test Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all tests
  python run_test_suite.py --all
  
  # Run specific category
  python run_test_suite.py --category unit
  
  # Run quick tests for development
  python run_test_suite.py --quick
  
  # Run with coverage report
  python run_test_suite.py --category unit --coverage
  
  # List available test categories
  python run_test_suite.py --list-categories
        """
    )
    
    parser.add_argument(
        "--category", "-c",
        choices=list(TestRunner(Path(__file__).parent.parent).test_categories.keys()),
        help="Run specific test category"
    )
    
    parser.add_argument(
        "--all", "-a",
        action="store_true",
        help="Run all test categories"
    )
    
    parser.add_argument(
        "--quick", "-q",
        action="store_true", 
        help="Run quick development tests"
    )
    
    parser.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        help="Exclude specific categories from --all"
    )
    
    parser.add_argument(
        "--list-categories",
        action="store_true",
        help="List available test categories"
    )
    
    parser.add_argument(
        "--check-env",
        action="store_true",
        help="Check test environment and dependencies"
    )
    
    parser.add_argument(
        "--install-deps",
        action="store_true",
        help="Install test dependencies"
    )
    
    parser.add_argument(
        "--coverage",
        action="store_true",
        help="Generate coverage report"
    )
    
    parser.add_argument(
        "--junit-xml",
        type=str,
        help="Save results in JUnit XML format"
    )
    
    parser.add_argument(
        "--timeout",
        type=int,
        default=3600,  # 1 hour
        help="Timeout for test execution in seconds"
    )
    
    args = parser.parse_args()
    
    # Initialize test runner
    project_root = Path(__file__).parent.parent
    runner = TestRunner(project_root)
    
    # Handle list categories
    if args.list_categories:
        print("Available test categories:")
        print("-" * 40)
        for category, config in runner.test_categories.items():
            print(f"{category:15} - {config.get('description', '')}")
        return 0
    
    # Handle environment check
    if args.check_env:
        success = runner.check_environment()
        return 0 if success else 1
    
    # Handle dependency installation
    if args.install_deps:
        success = runner.install_dependencies()
        return 0 if success else 1
    
    # Handle quick test
    if args.quick:
        success = runner.quick_test()
        return 0 if success else 1
    
    # Handle specific category
    if args.category:
        extra_args = []
        if args.coverage:
            extra_args.extend(["--cov-report=html", "--cov-report=term"])
        if args.junit_xml:
            extra_args.extend(["--junitxml", args.junit_xml])
        
        exit_code, duration, output = runner.run_test_category(args.category, extra_args)
        return exit_code
    
    # Handle run all
    if args.all:
        results = runner.run_all_categories(exclude=args.exclude)
        success = runner.generate_summary_report(results)
        return 0 if success else 1
    
    # Default: show help
    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())