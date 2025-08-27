# Test Organization Improvements for Better Parallelization

## Current Issues and Solutions

### 1. **Mixed Test Complexity**

**Problem**: Heavy training tests mixed with light unit tests in same files reduces parallel efficiency.

**Solution**: Reorganize tests into complexity-based directories:

```
tests/
├── unit/           # Fast unit tests (< 1s each)
│   ├── test_utils.py
│   ├── test_fft.py
│   ├── test_source.py
│   └── test_mrc.py
├── integration/    # Medium complexity (1-10s each)
│   ├── test_parse.py
│   ├── test_relion.py
│   └── test_writestar.py
├── training/       # Heavy training tests (> 10s each)
│   ├── test_reconstruct_fixed.py
│   ├── test_reconstruct_tilt.py
│   └── test_integration.py
└── data/           # Test data files
```

### 2. **Suboptimal Parametrization**

**Problem**: Large parametrized test matrices create many similar tests.

**Solution**: Use focused parametrization:

```python
# Instead of testing all combinations:
@pytest.mark.parametrize("particles", ["toy.mrcs", "toy.txt", "toy.star"])
@pytest.mark.parametrize("indices", [None, "first-100", "random-100"])
@pytest.mark.parametrize("ctf", [None, "CTF-Test"])

# Use focused combinations:
@pytest.mark.parametrize("particles,indices,ctf", [
    ("toy.mrcs", None, None),           # Basic case
    ("toy.mrcs", "first-100", "CTF-Test"), # Complex case
    ("toy.star", "random-100", None),   # Edge case
])
```

### 3. **Fixture Scope Optimization**

**Problem**: Function-scoped fixtures cause expensive repeated setup.

**Solution**: Use appropriate fixture scopes:

```python
@pytest.fixture(scope="session")  # For expensive one-time setup
def trained_model():
    # Train model once for entire session
    pass

@pytest.fixture(scope="class")    # For test class setup
def test_data():
    # Load data once per test class
    pass

@pytest.fixture(scope="function") # Only for test-specific setup
def temp_output_dir():
    # Create unique temp dir per test
    pass
```

## Optimized Test Execution Strategies

### Strategy 1: Stratified Execution (Recommended)

Run tests in phases based on complexity:

1. **Phase 1**: Fast unit tests with high parallelization (`-n4 --dist=each`)
2. **Phase 2**: Medium tests with moderate parallelization (`-n3 --dist=loadscope`) 
3. **Phase 3**: Heavy tests with conservative parallelization (`-n2 --dist=loadgroup`)

### Strategy 2: Resource-Aware Execution

Use custom test markers and execution:

```bash
# Run fast tests first
pytest -m "not slow" -n4 --dist=each

# Run slow tests separately  
pytest -m "slow" -n2 --dist=loadgroup
```

### Strategy 3: Test Splitting by Module

```bash
# Parallel execution of different test categories
pytest tests/unit/ -n4 --dist=each &
pytest tests/integration/ -n2 --dist=loadscope &
pytest tests/training/ -n1  # Sequential for heavy tests
wait
```

## Specific Improvements for Key Test Files

### test_integration.py
- **Current**: 13 test methods with expensive training
- **Improvement**: Split into `test_integration_light.py` and `test_integration_heavy.py`
- **Optimization**: Use session-scoped fixtures for trained models

### test_reconstruct_fixed.py  
- **Current**: 17 test methods with neural network training
- **Improvement**: Cache training results across similar parameter combinations
- **Optimization**: Use `@pytest.mark.slow` for resource management

### test_writestar.py
- **Current**: 4 test methods with heavy parametrization (48 combinations)
- **Improvement**: Reduce to focused parameter combinations (12 combinations)
- **Optimization**: Use class-scoped fixtures for common setup

## Performance Monitoring

Add test duration tracking:

```python
# In conftest.py
def pytest_runtest_teardown(item, nextitem):
    """Track test durations for optimization."""
    if hasattr(item, '_pytest_timing_start'):
        duration = time.time() - item._pytest_timing_start
        if duration > 10:  # Log slow tests
            print(f"SLOW TEST: {item.nodeid} took {duration:.2f}s")
```

## Resource Management

Prevent resource conflicts in parallel execution:

```python
# Use locks for heavy computational tests
@pytest.fixture(scope="session")
def computation_lock():
    import threading
    return threading.Lock()

def test_heavy_computation(computation_lock):
    with computation_lock:
        # Ensure only one heavy test runs at a time
        run_expensive_computation()
```

## Expected Performance Improvements

With these optimizations:

1. **Parallelization**: 4 CPUs instead of 2 → ~50% speed improvement
2. **Better load balancing**: `--dist=each` for mixed workloads → ~20% improvement  
3. **Reduced setup overhead**: Optimized fixture scopes → ~30% improvement
4. **Stratified execution**: Separate fast/slow tests → ~25% improvement

**Total expected improvement**: 60-80% reduction in test execution time.