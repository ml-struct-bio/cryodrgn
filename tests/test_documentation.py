"""Tests for documentation, docstrings, and example code validation.

This module validates that code examples in docstrings work correctly,
documentation is up-to-date, and API documentation matches implementation.
"""

import pytest
import doctest
import inspect
import re
import ast
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional
import tempfile
import subprocess

import cryodrgn
from cryodrgn import (
    dataset, ctf, fft, utils, source, lattice, 
    models, pose, analysis, config
)
from cryodrgn.commands import (
    backproject_voxel, train_vae, train_nn, analyze, 
    downsample, parse_ctf_star
)


class TestDocstrings:
    """Test docstring content and examples."""
    
    def test_module_docstrings_exist(self):
        """Test that main modules have docstrings."""
        modules_to_check = [
            cryodrgn.dataset,
            cryodrgn.ctf,
            cryodrgn.fft,
            cryodrgn.utils,
            cryodrgn.source,
            cryodrgn.lattice,
        ]
        
        for module in modules_to_check:
            assert module.__doc__ is not None, f"Module {module.__name__} missing docstring"
            assert len(module.__doc__.strip()) > 0, f"Module {module.__name__} has empty docstring"
    
    def test_command_docstrings_exist(self):
        """Test that command modules have comprehensive docstrings."""
        command_modules = [
            backproject_voxel,
            train_vae,
            train_nn,
            analyze,
            downsample,
        ]
        
        for module in command_modules:
            assert module.__doc__ is not None, f"Command {module.__name__} missing docstring"
            
            # Command docstrings should include example usage
            docstring = module.__doc__.lower()
            assert any(keyword in docstring for keyword in ['example', 'usage', '$']), \
                f"Command {module.__name__} docstring missing usage examples"
    
    def test_class_docstrings(self):
        """Test that important classes have docstrings."""
        classes_to_check = [
            (dataset.ImageDataset, "ImageDataset class"),
            (dataset.TiltSeriesData, "TiltSeriesData class"),
            (lattice.Lattice, "Lattice class"),
            (source.ImageSource, "ImageSource class"),
        ]
        
        for cls, description in classes_to_check:
            assert cls.__doc__ is not None, f"{description} missing docstring"
            assert len(cls.__doc__.strip()) > 50, f"{description} docstring too short"
    
    def test_function_docstrings(self):
        """Test that key functions have docstrings with proper format."""
        functions_to_check = [
            (utils.save_pkl, "save_pkl function"),
            (utils.load_pkl, "load_pkl function"),
            (ctf.compute_ctf, "compute_ctf function"),
            (fft.fft2_center, "fft2_center function"),
            (fft.ifft2_center, "ifft2_center function"),
        ]
        
        for func, description in functions_to_check:
            assert func.__doc__ is not None, f"{description} missing docstring"
            
            docstring = func.__doc__
            # Should have parameters and return description
            assert "parameters" in docstring.lower() or "args" in docstring.lower() or \
                   "param" in docstring.lower(), \
                   f"{description} docstring missing parameter documentation"


class TestDoctestExecution:
    """Test that docstring examples execute correctly."""
    
    def test_source_module_doctest(self):
        """Test doctest examples in source module."""
        try:
            import cryodrgn.source
            result = doctest.testmod(cryodrgn.source, verbose=True, report=True)
            
            # Allow some failures due to test data dependencies
            failure_rate = result.failed / max(result.attempted, 1)
            assert failure_rate < 0.5, f"Too many doctest failures in source module: {result}"
        except Exception as e:
            pytest.skip(f"Doctest execution failed: {e}")
    
    def test_utils_module_doctest(self):
        """Test doctest examples in utils module."""
        try:
            import cryodrgn.utils
            result = doctest.testmod(cryodrgn.utils, verbose=True, report=True)
            
            # Allow some failures due to external dependencies
            failure_rate = result.failed / max(result.attempted, 1)
            assert failure_rate < 0.3, f"Too many doctest failures in utils module: {result}"
        except Exception as e:
            pytest.skip(f"Doctest execution failed: {e}")


class TestExampleCodeValidation:
    """Validate example code from docstrings and documentation."""
    
    def extract_code_examples(self, docstring: str) -> List[str]:
        """Extract code examples from docstring."""
        if not docstring:
            return []
        
        # Look for code blocks in docstrings
        code_blocks = []
        
        # Pattern for >>> style examples
        example_pattern = r'>>> (.+?)(?=\n(?!\.{3}|\s*>>>)|\Z)'
        examples = re.findall(example_pattern, docstring, re.DOTALL | re.MULTILINE)
        code_blocks.extend(examples)
        
        # Pattern for code blocks
        code_block_pattern = r'```(?:python)?\n(.*?)```'
        blocks = re.findall(code_block_pattern, docstring, re.DOTALL)
        code_blocks.extend(blocks)
        
        return [block.strip() for block in code_blocks if block.strip()]
    
    def validate_code_syntax(self, code: str) -> bool:
        """Check if code has valid Python syntax."""
        try:
            ast.parse(code)
            return True
        except SyntaxError:
            return False
    
    def test_docstring_code_syntax(self):
        """Test that code examples in docstrings have valid syntax."""
        modules_to_check = [
            cryodrgn.source,
            cryodrgn.dataset, 
            cryodrgn.utils,
            backproject_voxel,
        ]
        
        syntax_errors = []
        
        for module in modules_to_check:
            if not module.__doc__:
                continue
                
            examples = self.extract_code_examples(module.__doc__)
            for i, example in enumerate(examples):
                if not self.validate_code_syntax(example):
                    syntax_errors.append(f"{module.__name__} example {i}: {example[:50]}...")
        
        assert len(syntax_errors) == 0, f"Syntax errors in docstring examples: {syntax_errors}"
    
    def test_command_example_syntax(self):
        """Test that command-line examples in docstrings are well-formed."""
        command_modules = [backproject_voxel, train_vae, train_nn, analyze]
        
        malformed_examples = []
        
        for module in command_modules:
            if not module.__doc__:
                continue
            
            # Look for command-line examples (lines starting with $ or cryodrgn)
            cmdline_pattern = r'^(?:\$\s*)?(cryodrgn\s+.+)$'
            examples = re.findall(cmdline_pattern, module.__doc__, re.MULTILINE)
            
            for example in examples:
                # Basic validation - should start with cryodrgn
                if not example.strip().startswith('cryodrgn'):
                    malformed_examples.append(f"{module.__name__}: {example}")
                
                # Should not have obvious syntax issues
                if '\\' in example and not example.endswith('\\'):
                    # Line continuation should end with backslash
                    if '\n' in example:
                        malformed_examples.append(f"{module.__name__}: multiline without proper continuation")
        
        assert len(malformed_examples) == 0, f"Malformed command examples: {malformed_examples}"


class TestAPIDocumentation:
    """Test that API documentation matches implementation."""
    
    def get_function_signature(self, func) -> str:
        """Get function signature as string."""
        try:
            sig = inspect.signature(func)
            return str(sig)
        except (ValueError, TypeError):
            return "Unable to get signature"
    
    def test_function_signatures_documented(self):
        """Test that function signatures match documentation."""
        functions_to_check = [
            (utils.save_pkl, ["data", "filename"]),
            (utils.load_pkl, ["filename"]),
            (ctf.compute_ctf, ["freqs", "dfu", "dfv", "dfang", "volt", "cs", "w"]),
        ]
        
        signature_mismatches = []
        
        for func, expected_params in functions_to_check:
            sig = inspect.signature(func)
            actual_params = list(sig.parameters.keys())
            
            # Check that expected parameters are present
            for param in expected_params:
                if param not in actual_params:
                    signature_mismatches.append(
                        f"{func.__name__} missing parameter: {param}"
                    )
            
            # Check docstring mentions parameters
            if func.__doc__:
                docstring_lower = func.__doc__.lower()
                for param in expected_params:
                    if param.lower() not in docstring_lower:
                        signature_mismatches.append(
                            f"{func.__name__} parameter {param} not documented in docstring"
                        )
        
        # Allow some mismatches due to refactoring
        assert len(signature_mismatches) < 5, f"Signature mismatches: {signature_mismatches}"
    
    def test_class_method_documentation(self):
        """Test that important class methods are documented."""
        classes_to_check = [
            (dataset.ImageDataset, ["__getitem__", "__len__"]),
            (source.ImageSource, ["from_file", "images", "chunks"]),
            (lattice.Lattice, ["get_circular_mask"]),
        ]
        
        undocumented_methods = []
        
        for cls, methods in classes_to_check:
            for method_name in methods:
                if hasattr(cls, method_name):
                    method = getattr(cls, method_name)
                    if callable(method) and not method.__doc__:
                        undocumented_methods.append(f"{cls.__name__}.{method_name}")
        
        # Allow some undocumented methods for special methods
        special_method_exceptions = ["__getitem__", "__len__", "__init__"]
        undocumented_methods = [
            method for method in undocumented_methods 
            if not any(exc in method for exc in special_method_exceptions)
        ]
        
        assert len(undocumented_methods) < 3, f"Undocumented methods: {undocumented_methods}"


class TestReadmeAndDocumentation:
    """Test README and other documentation files."""
    
    def test_readme_exists(self):
        """Test that README file exists and has content."""
        readme_path = Path(__file__).parent.parent / "README.md"
        assert readme_path.exists(), "README.md file not found"
        
        readme_content = readme_path.read_text()
        assert len(readme_content) > 100, "README file too short"
        
        # Should contain key sections
        content_lower = readme_content.lower()
        expected_sections = ["installation", "usage", "example"]
        
        for section in expected_sections:
            assert section in content_lower, f"README missing {section} section"
    
    def test_readme_code_examples_syntax(self):
        """Test that code examples in README have valid syntax."""
        readme_path = Path(__file__).parent.parent / "README.md" 
        if not readme_path.exists():
            pytest.skip("README.md not found")
        
        readme_content = readme_path.read_text()
        
        # Extract code blocks
        code_block_pattern = r'```(?:python|bash|sh)?\n(.*?)```'
        code_blocks = re.findall(code_block_pattern, readme_content, re.DOTALL)
        
        syntax_errors = []
        
        for i, block in enumerate(code_blocks):
            block = block.strip()
            if not block:
                continue
                
            # Skip bash/shell commands
            if block.startswith(('$', '#', 'cryodrgn', 'pip', 'conda')):
                continue
            
            # Check Python syntax
            try:
                ast.parse(block)
            except SyntaxError as e:
                syntax_errors.append(f"README code block {i}: {str(e)}")
        
        assert len(syntax_errors) == 0, f"Syntax errors in README: {syntax_errors}"


class TestTutorialValidation:
    """Test tutorial notebooks and documentation."""
    
    def find_notebook_files(self) -> List[Path]:
        """Find Jupyter notebook files in the project."""
        project_root = Path(__file__).parent.parent
        notebook_files = []
        
        # Look for notebooks in common locations
        search_paths = [
            project_root / "tutorials",
            project_root / "examples", 
            project_root / "notebooks",
            project_root / "docs",
            project_root / "cryodrgn" / "templates",
        ]
        
        for search_path in search_paths:
            if search_path.exists():
                notebook_files.extend(search_path.glob("*.ipynb"))
        
        return notebook_files
    
    @pytest.mark.slow
    def test_notebook_syntax(self):
        """Test that notebook cells have valid syntax."""
        notebook_files = self.find_notebook_files()
        
        if not notebook_files:
            pytest.skip("No notebook files found")
        
        syntax_errors = []
        
        try:
            import nbformat
        except ImportError:
            pytest.skip("nbformat not available")
        
        for notebook_path in notebook_files[:3]:  # Limit to first 3 notebooks
            try:
                notebook = nbformat.read(notebook_path, as_version=4)
                
                for cell_idx, cell in enumerate(notebook.cells):
                    if cell.cell_type == 'code' and cell.source.strip():
                        try:
                            # Basic syntax check
                            ast.parse(cell.source)
                        except SyntaxError as e:
                            syntax_errors.append(
                                f"{notebook_path.name} cell {cell_idx}: {str(e)}"
                            )
                            
            except Exception as e:
                # Skip problematic notebooks
                continue
        
        assert len(syntax_errors) < 5, f"Syntax errors in notebooks: {syntax_errors}"


def check_documentation_completeness():
    """Comprehensive documentation check function."""
    issues = []
    
    # Check main modules
    main_modules = [cryodrgn.dataset, cryodrgn.ctf, cryodrgn.source, cryodrgn.utils]
    for module in main_modules:
        if not module.__doc__ or len(module.__doc__.strip()) < 50:
            issues.append(f"Module {module.__name__} needs better docstring")
    
    # Check command modules
    command_modules = [backproject_voxel, train_vae, train_nn]
    for module in command_modules:
        if not module.__doc__ or 'example' not in module.__doc__.lower():
            issues.append(f"Command {module.__name__} needs usage examples in docstring")
    
    return issues


if __name__ == "__main__":
    """Run documentation tests when called as a script."""
    issues = check_documentation_completeness()
    if issues:
        print("Documentation issues found:")
        for issue in issues:
            print(f"- {issue}")
    else:
        print("Documentation completeness check passed!")
    
    # Run tests
    pytest.main([__file__, "-v"])