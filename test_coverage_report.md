# CryoDRGN Test Coverage Analysis Report
==================================================

## Test Suite Statistics
📁 Total test files: 35
🧪 Total test functions: 163
📊 Parametrized tests: 222
🏗️ Class-based tests: 23
🚀 Training tests: 9
💾 File I/O tests: 24

## Test Organization Metrics
📈 Average tests per file: 4.7
🎯 Parametrization ratio: 136.2%

## Potential Coverage Gaps
⚠️ Modules without dedicated tests (44 total):
   - _version
   - abinit_het
   - abinit_homo
   - analysis
   - analyze
   - analyze_convergence
   - analyze_landscape
   - analyze_landscape_full
   - backproject_voxel
   - beta_schedule
   ... and 34 more

## Improvement Suggestions
- 🔍 Add tests for untested modules: _version, abinit_het, abinit_homo, analysis, analyze
- ⚡ Consider optimizing 9 training tests with session-scoped fixtures
- 📊 High number of parametrized tests - consider reducing combinations for faster execution

## Parallelization Opportunities
- Fast unit tests: FFT, utils, source operations
- Medium tests: Parsing, file I/O, data processing
- Heavy tests: Training, reconstruction, integration

## Recommended Test Execution Strategy
```bash
# Phase 1: Fast tests (high parallelization)
pytest tests/test_utils.py tests/test_fft.py tests/test_source.py -n4 --dist=each

# Phase 2: Medium tests (moderate parallelization)
pytest tests/test_parse.py tests/test_relion.py tests/test_writestar.py -n3 --dist=loadscope

# Phase 3: Heavy tests (conservative parallelization)
pytest tests/test_integration.py tests/test_reconstruct_*.py -n2 --dist=loadgroup
```