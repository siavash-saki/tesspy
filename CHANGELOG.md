# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-02-22

### Added
- `pyproject.toml` replaces `setup.py` (Hatchling build backend)
- `tesspy/_version.py` — single source of truth for the package version
- `tesspy/_constants.py` — `OSM_PRIMARY_FEATURES` and `OSM_HIGHWAY_TYPES` constants
  shared across the package (eliminates duplication across three files)
- `tesspy/_validators.py` — input validation helpers extracted from `tessellation.py`
- `tesspy/methods/` sub-package with one module per tessellation algorithm:
  - `squares.py` — `get_squares_polyfill`, `get_adaptive_squares`, `count_poi`
  - `hexagons.py` — `get_h3_hexagons`
  - `voronoi.py` — `voronoi_polygons`
  - `city_blocks.py` — `split_linestring`, `create_blocks`, `get_rest_polygon`, `explode`
  - `_clustering.py` — `get_hierarchical_clustering_parameter`
- `tesspy/data/` sub-package for data retrieval:
  - `poi.py` — `POIdata` class
  - `roads.py` — `RoadData` class
  - `_overpass.py` — internal coordinate helpers
  - `_geo.py` — `get_city_polygon`, `count_poi_per_tile`
- Type annotations on all public API surfaces
- `.pre-commit-config.yaml` with ruff (lint + format) and common hooks
- `tests/conftest.py` with session-scoped fixtures and `call_with_osm_retry` helper
- `tests/unit/` — fast tests requiring no network access
- `tests/integration/` — OSM API tests marked with `@pytest.mark.integration`
- `CHANGELOG.md` (this file)

### Changed
- Minimum Python version raised from 3.7 to 3.9 (3.7 and 3.8 are EOL)
- `tessellation.py` reduced from ~785 lines to ~250 lines (methods delegate to `tesspy/methods/`)
- CI (`tests_package.yml`) updated:
  - Replaced conda/Mambaforge with plain `pip install -e ".[dev]"`
  - Unit tests run on every push/PR across Python 3.9–3.12
  - Integration tests run only on pushes to `main` (saves OSM API quota)
  - Replaced `psf/black@stable` with ruff (lint + format)
  - Added mypy type checking step
  - Updated all action versions to v4/v5
- `python-publish.yml` switched from username/password secrets to OIDC trusted publishing
  (configure trusted publishing in PyPI project settings to activate)

### Fixed
- Version inconsistency: `__init__.py` reported `0.0.1` while `setup.py` reported `0.1.2`
- Relative path bug in `test_tessellation_functions.py` (`"tests/..."` → `Path(__file__).parent`)
- Wildcard imports (`from tessellation_functions import *`) replaced with explicit imports

### Deprecated
- `tesspy/tessellation_functions.py` as a direct import source — now a shim;
  use `tesspy.methods.*` in new code
- `tesspy/poi_data.py` as a direct import source — now a shim;
  use `tesspy.data.*` in new code

---

## [0.1.2] — 2022

### Changed
- Minor bug fixes and improvements

## [0.1.0] — 2022

### Added
- Initial public release with five tessellation methods:
  squares, hexagons, adaptive squares, Voronoi, city blocks
- Published in Journal of Open Source Software (JOSS), DOI: 10.21105/joss.04620
