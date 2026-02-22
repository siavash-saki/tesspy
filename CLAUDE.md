# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

tesspy is a Python library for geographical tessellation using OpenStreetMap data. It provides five tessellation methods (squares, adaptive squares, hexagons, Voronoi, city blocks) over arbitrary study areas defined by city name or GeoDataFrame. Python 3.11+ required.

## Common Commands

```bash
# Install in dev mode (with example notebook deps)
uv pip install -e ".[dev,examples]"

# Unit tests (no network calls)
pytest tests/unit/ -v

# Single test file
pytest tests/unit/test_squares.py -v

# Integration tests (hits OSM API, use sparingly)
pytest tests/integration/ -v -m integration --timeout=600

# Coverage
pytest tests/unit/ --cov=tesspy --cov-report=term-missing

# Lint and format
ruff check tesspy/
ruff format tesspy/

# Type checking
mypy tesspy/ --ignore-missing-imports

# Pre-commit hooks
pre-commit run --all-files
```

## Architecture

Three-layer design:

**Data Layer** (`tesspy/data/`) — All external I/O (OSM Overpass API, osmnx, geocoding). `POIdata` for POI queries, `RoadData` for road networks, `get_city_polygon()` for geocoding.

**Algorithm Layer** (`tesspy/methods/`) — Pure computational functions, no I/O. Five modules: `squares.py`, `hexagons.py`, `voronoi.py`, `city_blocks.py`, `_clustering.py`.

**Orchestration Layer** (`tesspy/tessellation.py`) — `Tessellation` class is the main user-facing interface. Accepts city name (string) or GeoDataFrame, exposes methods like `.squares()`, `.hexagons()`, `.voronoi()`, `.city_blocks()`, `.adaptive_squares()`.

**Public API** (exported from `__init__.py`): `Tessellation`, `count_poi_per_tile`, `configure_logging`, `__version__`.

## Key Conventions

- **Logging**: Package logger is `logging.getLogger("tesspy")`. Silent by default (NullHandler). Users enable via `configure_logging()`. Use `log_progress()` helper for verbose output — never use `print()`.
- **Constants**: `_constants.py` is single source of truth for `DEFAULT_POI_CATEGORIES`, `OSM_PRIMARY_FEATURES`, `OSM_HIGHWAY_TYPES`.
- **Validation**: Shared validators in `_validators.py` (`_check_input_geodataframe`, `_check_valid_geometry_gdf`).
- **Type hints**: All public APIs must have type annotations. mypy enforced in CI.
- **Style**: ruff with rules E, F, W, I, UP, B, C4. Line length 88.
- **Testing**: pytest with `@pytest.mark.integration` for OSM API tests, `@pytest.mark.slow` for long tests. Unit tests use fixtures from `tests/conftest.py` and `tests/cache/` (no network). Use `call_with_osm_retry()` from conftest for integration test resilience.

## Best Practices

- Always suggest commit messages after editing files.

## CI Pipeline

Runs on push to main and PRs to main: ruff lint/format check → mypy → unit tests (Python 3.11/3.12/3.13 matrix) → integration tests (main branch only, to preserve OSM API quota).

## Branch Strategy

`main` is stable release. `develop` is active development. Feature branches from `develop`.
