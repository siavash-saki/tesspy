# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-02-22

### Added
- `pyproject.toml` replaces `setup.py` (Hatchling build backend)
- `tesspy/_version.py` — single source of truth for the package version
- `tesspy/_constants.py` — `DEFAULT_POI_CATEGORIES`, `OSM_PRIMARY_FEATURES`, and
  `OSM_HIGHWAY_TYPES` constants shared across the package (eliminates duplication)
- `tesspy/_validators.py` — input validation helpers extracted from `tessellation.py`
- `tesspy/methods/` sub-package with one module per tessellation algorithm:
  - `squares.py` — `get_squares_polyfill`, `get_adaptive_squares`, `count_poi`
  - `hexagons.py` — `get_h3_hexagons`
  - `voronoi.py` — `voronoi_polygons`
  - `city_blocks.py` — `create_city_blocks` (new primary API),
    `merge_city_blocks` (new), plus legacy `split_linestring`, `create_blocks`,
    `get_rest_polygon`, `explode`
  - `_clustering.py` — `get_hierarchical_clustering_parameter`
- `tesspy/data/` sub-package for data retrieval:
  - `poi.py` — `POIdata` class
  - `roads.py` — `RoadData` class
  - `_overpass.py` — coordinate rounding helpers (`geom_ceil`, `geom_floor`)
  - `_geo.py` — `get_city_polygon`, `count_poi_per_tile`
- `tesspy/_logging.py` — `configure_logging()` and `log_progress()` helpers;
  package auto-configures at INFO level on import
- Type annotations on all public API surfaces
- `.pre-commit-config.yaml` with ruff (lint + format) and common hooks
- `tests/conftest.py` with session-scoped fixtures and `call_with_osm_retry` helper
- `tests/unit/` — fast tests requiring no network access
- `tests/integration/` — OSM API tests marked with `@pytest.mark.integration`
- `uv.lock` — lock file for reproducible installs with uv
- `CHANGELOG.md` (this file)
- Docs: Sphinx modernized with MyST parser, Furo theme, and deduplicated example
  notebooks

### Changed
- Minimum Python version raised from **3.7 to 3.11** (3.7–3.10 are EOL or unsupported
  by current dependency minimums)
- `tessellation.py` reduced from ~784 lines to ~530 lines (methods delegate to
  `tesspy/methods/`)
- `POIdata` now uses osmnx for data retrieval instead of raw `requests`/JSON Overpass
  calls; benefits from osmnx's automatic retry, rate-limit handling, and caching
- `city_blocks`: boundary-inclusive polygonization with adjacency-constrained merging
  via new `create_city_blocks` / `merge_city_blocks` API
- Codebase modernized for latest pandas (≥2.0), geopandas (≥1.0), and scikit-learn
  (≥1.3) APIs — no deprecated method calls remain
- All `print()` statements replaced with structured Python logging
- Default Overpass API timeout increased from 60 s to 300 s
- `hdbscan` package dependency **replaced** by `scikit-learn>=1.3`, which ships its
  own `HDBSCAN` implementation; the `cluster_algo='hdbscan'` option in `.voronoi()`
  now uses `sklearn.cluster.HDBSCAN`
- `matplotlib` added as an explicit runtime dependency
- All runtime dependencies now carry explicit minimum versions:
  `geopandas>=1.0`, `scipy>=1.11`, `h3>=4.0`, `osmnx>=1.9`, `mercantile>=1.2`,
  `scikit-learn>=1.3`, `requests>=2.28`, `numpy>=1.26`, `pandas>=2.0`,
  `matplotlib>=3.10.8`
- CI (`tests_package.yml`) updated:
  - Triggers on push/PR to both `main` and `develop` branches
  - Replaced conda/Mambaforge with plain `pip install -e ".[dev]"`
  - Lint job uses `uvx ruff` (no Python install needed)
  - Unit tests run on every push/PR on **Python 3.11, 3.12, 3.13** (ubuntu only)
  - Integration tests run only on pushes to `main` (saves OSM API quota)
  - Replaced `psf/black@stable` with ruff (lint + format)
  - Added mypy type-checking job
  - Added Codecov coverage upload
  - Updated all action versions to v4/v5
- `python-publish.yml` switched from Twine + username/password secrets to OIDC trusted
  publishing via `pypa/gh-action-pypi-publish` (configure trusted publishing in PyPI
  project settings to activate)
- Installation docs updated: uv is now the recommended install method

### Performance
- Squares: batch DataFrame operations instead of per-row concat
- Hexagons: unified Polygon/MultiPolygon handling into a single GeoDataFrame path
- POI: optimized category matching and node geometry creation
- City blocks: vectorized linestring splitting and geometry validation
- Clustering: binary search for hierarchical clustering distance threshold
- `count_poi_per_tile`: replaced `pivot_table` with `pd.crosstab`
- `Tessellation`: replaced per-tile loops with vectorized/grouped operations

### Fixed
- Version inconsistency: `__init__.py` reported `0.0.1` while `setup.py` reported
  `0.1.2`
- Relative path bug in `test_tessellation_functions.py` (`"tests/..."` →
  `Path(__file__).parent`)
- Wildcard imports (`from tessellation_functions import *`) replaced with explicit
  imports
- h3 API migration: `hexagons` module updated to h3 v4 API
- CRS assignment: `set_crs()` used instead of direct attribute assignment
- Index duplication after GeoDataFrame joins (`.reset_index(drop=True)`)
- mypy type errors resolved across all modules
- ruff lint/format issues resolved
- Hatchling `pyproject.toml` typo corrected
- Notebook and library compatibility issues resolved
- `assert` removed from runtime code paths (ruff `B011`)
- Deprecated pandas method in test suite corrected

### Deprecated
- `tesspy/tessellation_functions.py` as a direct import source — now a shim that
  emits `DeprecationWarning`; **will be removed in v0.3.0**.
  Use `tesspy.methods.*` in new code.
- `tesspy/poi_data.py` as a direct import source — now a shim that emits
  `DeprecationWarning`; **will be removed in v0.3.0**.
  Use `tesspy.data.*` in new code.

### Removed
- `setup.py` — replaced by `pyproject.toml`
- `hdbscan` package from runtime dependencies — covered by `scikit-learn>=1.3`

---

## [0.1.2] — 2022

### Changed
- Minor bug fixes and improvements

## [0.1.0] — 2022

### Added
- Initial public release with five tessellation methods:
  squares, hexagons, adaptive squares, Voronoi, city blocks
- Published in Journal of Open Source Software (JOSS), DOI: 10.21105/joss.04620
