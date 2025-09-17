# cryoDRGN Testing and Development Makefile

.PHONY: help install test test-all test-unit test-integration test-performance test-quick test-coverage clean lint format check-env

# Default target
help: ## Show this help message
	@echo "cryoDRGN Development Commands"
	@echo "============================="
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install cryoDRGN with development dependencies
	python -m pip install -e .[dev]

install-test-deps: ## Install additional test dependencies
	python tests/run_test_suite.py --install-deps

check-env: ## Check testing environment and dependencies
	python tests/run_test_suite.py --check-env

test: ## Run default test suite (unit tests)
	python tests/run_test_suite.py --category unit

test-all: ## Run all test categories
	python tests/run_test_suite.py --all

test-unit: ## Run unit tests only
	python tests/run_test_suite.py --category unit

test-integration: ## Run integration tests
	python tests/run_test_suite.py --category integration

test-performance: ## Run performance benchmark tests
	python tests/run_test_suite.py --category performance

test-error-handling: ## Run error handling and edge case tests
	python tests/run_test_suite.py --category error_handling

test-property-based: ## Run property-based tests with Hypothesis
	python tests/run_test_suite.py --category property_based

test-mocking: ## Run tests with mocked dependencies
	python tests/run_test_suite.py --category mocking

test-documentation: ## Run documentation and docstring tests
	python tests/run_test_suite.py --category documentation

test-quick: ## Run quick development tests
	python tests/run_test_suite.py --quick

test-coverage: ## Run tests with coverage report
	python tests/run_test_suite.py --category unit --coverage
	@echo "Coverage report available in htmlcov/index.html"

test-slow: ## Run slow/long-running tests
	python tests/run_test_suite.py --category slow

test-gpu: ## Run GPU-specific tests (requires GPU)
	python tests/run_test_suite.py --category gpu

test-cmdline: ## Run command-line interface tests
	python tests/run_test_suite.py --category cmdline

list-test-categories: ## List all available test categories
	python tests/run_test_suite.py --list-categories

# Legacy pytest commands for compatibility
pytest: ## Run pytest with basic configuration
	python -m pytest tests/ -v

pytest-coverage: ## Run pytest with coverage
	python -m pytest tests/ -v --cov=cryodrgn --cov-report=html --cov-report=term

pytest-quick: ## Run pytest with fail-fast for development
	python -m pytest tests/ -v --maxfail=3 -x

# Code quality
lint: ## Run linting checks
	pre-commit run --all-files

format: ## Format code with black and other formatters  
	pre-commit run --all-files black
	pre-commit run --all-files isort

check: ## Run all code quality checks
	pre-commit run --all-files
	python tests/run_test_suite.py --category documentation

# Cleanup
clean: ## Clean up temporary files and test artifacts
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	rm -rf build/
	rm -rf dist/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf test_results/
	rm -rf .pytest_cache/
	rm -rf .hypothesis/

clean-test-results: ## Clean up test result files
	rm -rf test_results/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf coverage.xml
	rm -rf .pytest_cache/

# Development helpers
dev-setup: install check-env ## Set up development environment
	@echo "Development environment ready!"
	@echo "Run 'make test-quick' to verify installation"

benchmark: ## Run performance benchmarks and save results
	python tests/run_test_suite.py --category performance
	@echo "Benchmark results saved in test_results/"

# CI simulation
ci-test: ## Simulate CI testing locally
	python tests/run_test_suite.py --all --exclude slow performance
	@echo "CI simulation complete"

# Documentation
docs-test: ## Test documentation examples
	python tests/run_test_suite.py --category documentation

# Docker support (if Dockerfile exists)
docker-test: ## Run tests in Docker container
	@if [ -f Dockerfile ]; then \
		docker build -t cryodrgn-test . && \
		docker run --rm cryodrgn-test make test-all; \
	else \
		echo "No Dockerfile found"; \
	fi

# Performance monitoring
profile-tests: ## Profile test performance
	python -m cProfile -o test_profile.prof tests/run_test_suite.py --category performance
	@echo "Profile saved to test_profile.prof"

# Git hooks
install-hooks: ## Install git pre-commit hooks
	pre-commit install

# Dependencies
update-deps: ## Update development dependencies
	python -m pip install --upgrade pip
	python -m pip install --upgrade -e .[dev]

# Help for specific test categories
help-testing: ## Show detailed testing help
	@echo ""
	@echo "Testing Categories:"
	@echo "=================="
	@python tests/run_test_suite.py --list-categories
	@echo ""
	@echo "Examples:"
	@echo "  make test-unit          # Fast unit tests"
	@echo "  make test-integration   # Integration tests"
	@echo "  make test-quick         # Quick dev tests"
	@echo "  make test-all           # All test categories"
	@echo "  make test-coverage      # Tests with coverage"
	@echo ""