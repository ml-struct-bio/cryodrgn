#!/usr/bin/env python3
"""
Test coverage analysis script for cryodrgn_beta.

This script analyzes the current test suite to identify potential gaps in coverage
and suggests improvements.
"""

import os
import ast
import re
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple


class CoverageAnalyzer:
    """Analyze test coverage patterns and identify gaps."""
    
    def __init__(self, src_dir: str = "cryodrgn", test_dir: str = "tests"):
        self.src_dir = Path(src_dir)
        self.test_dir = Path(test_dir)
        self.src_files = list(self.src_dir.rglob("*.py"))
        self.test_files = list(self.test_dir.rglob("test_*.py"))
        
    def find_source_functions(self) -> Dict[str, List[str]]:
        """Find all functions/classes in source code."""
        functions = defaultdict(list)
        
        for src_file in self.src_files:
            if "__pycache__" in str(src_file):
                continue
                
            try:
                with open(src_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                tree = ast.parse(content)
                
                for node in ast.walk(tree):
                    if isinstance(node, ast.FunctionDef):
                        if not node.name.startswith('_'):  # Skip private functions
                            functions[str(src_file.relative_to(self.src_dir))].append(
                                f"function:{node.name}"
                            )
                    elif isinstance(node, ast.ClassDef):
                        functions[str(src_file.relative_to(self.src_dir))].append(
                            f"class:{node.name}"
                        )
                        
            except Exception as e:
                print(f"Error parsing {src_file}: {e}")
                
        return functions
    
    def find_tested_items(self) -> Set[str]:
        """Find what items are being tested based on import patterns and test names."""
        tested_items = set()
        
        for test_file in self.test_files:
            try:
                with open(test_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # Find imports from cryodrgn
                import_matches = re.findall(r'from cryodrgn\.(\S+) import (.+)', content)
                for module, items in import_matches:
                    for item in items.split(','):
                        item = item.strip()
                        if item and item != '*':
                            tested_items.add(f"{module}.py:function:{item}")
                
                # Find test function names that might indicate what's being tested
                test_functions = re.findall(r'def (test_\w+)', content)
                for func in test_functions:
                    # Infer what might be tested from test name
                    tested_name = func.replace('test_', '').replace('_', '')
                    tested_items.add(f"inferred:{tested_name}")
                    
            except Exception as e:
                print(f"Error analyzing {test_file}: {e}")
                
        return tested_items
    
    def analyze_test_patterns(self) -> Dict[str, int]:
        """Analyze patterns in test organization."""
        patterns = defaultdict(int)
        
        for test_file in self.test_files:
            try:
                with open(test_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # Count test organization patterns
                patterns['total_test_files'] += 1
                patterns['total_test_functions'] += len(re.findall(r'def test_', content))
                patterns['parametrized_tests'] += len(re.findall(r'@pytest\.mark\.parametrize', content))
                patterns['class_based_tests'] += len(re.findall(r'class Test\w+', content))
                patterns['fixture_usage'] += len(re.findall(r'def \w+\([^)]*\w+\)', content))
                
                # Check for expensive operations
                if 'train' in content.lower() or 'epoch' in content.lower():
                    patterns['training_tests'] += 1
                if 'tmpdir' in content or 'tmp_path' in content:
                    patterns['file_io_tests'] += 1
                    
            except Exception as e:
                print(f"Error analyzing {test_file}: {e}")
                
        return patterns
    
    def identify_untested_modules(self) -> List[str]:
        """Identify source modules that may lack test coverage."""
        src_modules = {f.stem for f in self.src_files if f.stem != '__init__'}
        test_modules = {f.stem.replace('test_', '') for f in self.test_files}
        
        # Also check for modules referenced in test files
        referenced_modules = set()
        for test_file in self.test_files:
            try:
                with open(test_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                    
                # Find module references
                modules = re.findall(r'from cryodrgn\.(\w+)', content)
                referenced_modules.update(modules)
                modules = re.findall(r'import cryodrgn\.(\w+)', content)
                referenced_modules.update(modules)
                
            except Exception:
                continue
        
        untested = src_modules - test_modules - referenced_modules
        return sorted(list(untested))
    
    def suggest_improvements(self) -> List[str]:
        """Suggest specific improvements for test coverage and performance."""
        suggestions = []
        
        patterns = self.analyze_test_patterns()
        untested = self.identify_untested_modules()
        
        # Coverage suggestions
        if untested:
            suggestions.append(f"🔍 Add tests for untested modules: {', '.join(untested[:5])}")
            
        # Performance suggestions
        if patterns['training_tests'] > 3:
            suggestions.append(f"⚡ Consider optimizing {patterns['training_tests']} training tests with session-scoped fixtures")
            
        if patterns['parametrized_tests'] > 50:
            suggestions.append("📊 High number of parametrized tests - consider reducing combinations for faster execution")
            
        # Organization suggestions
        ratio = patterns['total_test_functions'] / max(patterns['total_test_files'], 1)
        if ratio > 20:
            suggestions.append("🗂️ Consider splitting large test files for better parallelization")
            
        if patterns['class_based_tests'] < patterns['total_test_files'] * 0.3:
            suggestions.append("🏗️ Consider using more class-based tests for better fixture management")
            
        return suggestions
    
    def generate_report(self) -> str:
        """Generate a comprehensive coverage analysis report."""
        patterns = self.analyze_test_patterns()
        untested = self.identify_untested_modules()
        suggestions = self.suggest_improvements()
        
        report = [
            "# CryoDRGN Test Coverage Analysis Report",
            "=" * 50,
            "",
            "## Test Suite Statistics",
            f"📁 Total test files: {patterns['total_test_files']}",
            f"🧪 Total test functions: {patterns['total_test_functions']}",
            f"📊 Parametrized tests: {patterns['parametrized_tests']}",
            f"🏗️ Class-based tests: {patterns['class_based_tests']}",
            f"🚀 Training tests: {patterns['training_tests']}",
            f"💾 File I/O tests: {patterns['file_io_tests']}",
            "",
            f"## Test Organization Metrics",
            f"📈 Average tests per file: {patterns['total_test_functions'] / max(patterns['total_test_files'], 1):.1f}",
            f"🎯 Parametrization ratio: {patterns['parametrized_tests'] / max(patterns['total_test_functions'], 1):.1%}",
            "",
            "## Potential Coverage Gaps",
        ]
        
        if untested:
            report.extend([
                f"⚠️ Modules without dedicated tests ({len(untested)} total):",
                *[f"   - {module}" for module in untested[:10]],
            ])
            if len(untested) > 10:
                report.append(f"   ... and {len(untested) - 10} more")
        else:
            report.append("✅ All major modules appear to have test coverage")
            
        report.extend([
            "",
            "## Improvement Suggestions",
            *[f"- {suggestion}" for suggestion in suggestions],
            "",
            "## Parallelization Opportunities",
            "- Fast unit tests: FFT, utils, source operations",
            "- Medium tests: Parsing, file I/O, data processing", 
            "- Heavy tests: Training, reconstruction, integration",
            "",
            "## Recommended Test Execution Strategy",
            "```bash",
            "# Phase 1: Fast tests (high parallelization)",
            "pytest tests/test_utils.py tests/test_fft.py tests/test_source.py -n4 --dist=each",
            "",
            "# Phase 2: Medium tests (moderate parallelization)",
            "pytest tests/test_parse.py tests/test_relion.py tests/test_writestar.py -n3 --dist=loadscope",
            "",
            "# Phase 3: Heavy tests (conservative parallelization)",
            "pytest tests/test_integration.py tests/test_reconstruct_*.py -n2 --dist=loadgroup",
            "```",
        ])
        
        return "\n".join(report)


def main():
    """Run the coverage analysis."""
    print("🔍 Analyzing cryoDRGN test coverage...")
    
    analyzer = CoverageAnalyzer()
    report = analyzer.generate_report()
    
    # Save report to file
    with open("test_coverage_report.md", "w") as f:
        f.write(report)
    
    print("✅ Analysis complete!")
    print("\n" + report)
    print(f"\n📄 Detailed report saved to: test_coverage_report.md")


if __name__ == "__main__":
    main()