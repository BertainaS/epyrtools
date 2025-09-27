# Makefile for EPyR Tools development and testing

.PHONY: help install install-dev test test-cov lint format type-check security clean build docs
.DEFAULT_GOAL := help

# Variables
PYTHON := python
PIP := pip
PYTEST := pytest
PACKAGE := epyr

help: ## Show this help message
	@echo "EPyR Tools - Development Commands"
	@echo "================================="
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install: ## Install package in development mode
	$(PIP) install -e .

install-dev: ## Install package with development dependencies
	$(PIP) install -e ".[dev,docs]"
	pre-commit install

# Testing targets
test: ## Run tests
	$(PYTEST) tests/ -v

test-cov: ## Run tests with coverage
	$(PYTEST) tests/ -v --cov=$(PACKAGE) --cov-report=html --cov-report=term-missing

test-fast: ## Run fast tests (exclude slow markers)
	$(PYTEST) tests/ -v -m "not slow"

test-integration: ## Run integration tests only
	$(PYTEST) tests/test_integration.py -v

test-core: ## Run core functionality tests for release validation
	$(PYTEST) tests/test_core_functionality.py -v

# Code quality targets
lint: ## Run linting (flake8)
	flake8 $(PACKAGE)/ tests/

format: ## Format code with black and isort
	black $(PACKAGE)/ tests/
	isort $(PACKAGE)/ tests/

format-check: ## Check code formatting without making changes
	black --check $(PACKAGE)/ tests/
	isort --check-only $(PACKAGE)/ tests/

type-check: ## Run type checking with mypy
	mypy $(PACKAGE)/

security: ## Run security checks with bandit
	bandit -r $(PACKAGE)/ -c pyproject.toml

# Quality checks (all at once)
quality: lint type-check security ## Run all quality checks

quality-fix: format lint ## Fix formatting and run linting

# Pre-commit
pre-commit: ## Run pre-commit hooks on all files
	pre-commit run --all-files

# Documentation targets  
docs: ## Build documentation
	cd docs && $(MAKE) html

docs-clean: ## Clean documentation build
	cd docs && $(MAKE) clean

docs-serve: docs ## Build and serve documentation locally
	cd docs/_build/html && $(PYTHON) -m http.server 8000

# Package building
build: clean ## Build package
	$(PYTHON) -m build

build-wheel: clean ## Build wheel only
	$(PYTHON) -m build --wheel

build-sdist: clean ## Build source distribution only
	$(PYTHON) -m build --sdist

# CLI testing
cli-test: ## Test CLI commands
	epyr --help
	epyr-info --help
	epyr-config show plotting

# Cleanup targets
clean: ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

clean-all: clean ## Clean everything including virtual environment
	rm -rf .tox/
	rm -rf .venv/

# Development environment
venv: ## Create virtual environment
	$(PYTHON) -m venv .venv
	@echo "Activate with: source .venv/bin/activate"

venv-install: venv ## Create venv and install package
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e ".[dev,docs]"
	.venv/bin/pre-commit install

# Performance and benchmarks
benchmark: ## Run performance benchmarks
	$(PYTEST) tests/ -v -m "benchmark"

profile: ## Run profiling tests
	$(PYTEST) tests/ -v --profile

# Release helpers
version: ## Show current version
	$(PYTHON) -c "import epyr; print(f'Version: {epyr.__version__}')"

check-version: ## Check if version is consistent across files
	@echo "Checking version consistency..."
	@VERSION_INIT=$$($(PYTHON) -c "import epyr; print(epyr.__version__)"); \
	VERSION_TOML=$$(grep '^version = ' pyproject.toml | cut -d'"' -f2); \
	if [ "$$VERSION_INIT" != "$$VERSION_TOML" ]; then \
		echo "Version mismatch: __init__.py=$$VERSION_INIT, pyproject.toml=$$VERSION_TOML"; \
		exit 1; \
	else \
		echo "Version $$VERSION_INIT is consistent"; \
	fi

# CI/CD simulation
ci: install-dev quality test-cov build ## Simulate CI pipeline locally

# Development workflow
dev-setup: venv-install ## Complete development setup
	@echo "Development environment ready!"
	@echo "Run 'source .venv/bin/activate' to activate virtual environment"

# Advanced testing
test-all: test-cov test-integration benchmark ## Run all tests including benchmarks

# Maintenance
update-deps: ## Update development dependencies
	$(PIP) install --upgrade pip
	$(PIP) install --upgrade -e ".[dev,docs]"
	pre-commit autoupdate

# Docker support (if needed)
docker-test: ## Run tests in Docker container
	docker build -t epyr-test -f Dockerfile.test .
	docker run --rm epyr-test

# Information targets
info: ## Show project information
	@echo "EPyR Tools - Project Information"
	@echo "==============================="
	@echo "Python version: $$($(PYTHON) --version)"
	@echo "Package location: $$($(PYTHON) -c 'import epyr; print(epyr.__file__)')"
	@echo "Git branch: $$(git branch --show-current 2>/dev/null || echo 'Not in git repo')"
	@echo "Git status: $$(git status --porcelain | wc -l) modified files"

deps-tree: ## Show dependency tree
	pipdeptree

# Performance monitoring
memory-test: ## Run memory usage tests
	$(PYTEST) tests/test_performance.py -v -k "memory"

benchmark-all: ## Run comprehensive benchmarks
	$(PYTEST) tests/ -v -m "benchmark" --benchmark-sort=mean