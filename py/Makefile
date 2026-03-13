# ============================================================
# spintronics Python bindings — Makefile
# COOLJAPAN OU (Team Kitasan)
# ============================================================

.PHONY: build build-debug install dev test clean wheel lint fmt check

# Default: build release wheel
build:
	maturin build --release

# Fast debug / editable install (no optimisations)
build-debug:
	maturin develop

# Editable install — release optimisations
dev: install

install:
	maturin develop --release

# Build a redistributable wheel into dist/
wheel:
	maturin build --release --out dist

# Run Rust tests with the python feature enabled
test:
	cd .. && cargo nextest run --features python --no-fail-fast

# Run Python-side tests (requires maturin develop first)
pytest:
	pytest tests/ -v

# Lint Python code
lint:
	flake8 python/spintronics/ tests/ --max-line-length 88
	isort --check-only --diff python/spintronics/ tests/

# Format Python code
fmt:
	black python/spintronics/ tests/
	isort python/spintronics/ tests/

# Type-check Python stubs
check:
	mypy python/spintronics/ --ignore-missing-imports

# Remove build artefacts
clean:
	rm -rf dist/ ../target/wheels/ __pycache__ .pytest_cache
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

# Show help
help:
	@echo "spintronics Python bindings — available targets:"
	@echo "  build        Build release wheel via maturin"
	@echo "  build-debug  Editable debug install (maturin develop)"
	@echo "  install/dev  Editable release install (maturin develop --release)"
	@echo "  wheel        Build redistributable .whl into dist/"
	@echo "  test         Run Rust tests with --features python"
	@echo "  pytest       Run Python test suite"
	@echo "  lint         Flake8 + isort checks"
	@echo "  fmt          Black + isort autoformat"
	@echo "  check        mypy type-check stubs"
	@echo "  clean        Remove dist/, wheels, pycache artefacts"
