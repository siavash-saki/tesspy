# tesspy v0.2.0 — Improvement Roadmap

This roadmap covers all remaining improvements for the v0.2.0 release, starting from the `develop` branch. Group 1 (architecture & modularization) is already complete — see [CHANGELOG.md](CHANGELOG.md) for details.

---

## Group 2 — Code Quality & Bug Fixes

Fix type-checking anti-patterns, unsafe attribute lookups, and warning suppression issues.

**Files:** `tesspy/tessellation.py`, `tesspy/methods/hexagons.py`, `tesspy/data/roads.py`, `tesspy/data/poi.py`, `tesspy/data/_geo.py`, `pyproject.toml`

- [x] Replace `type(area) == gpd.GeoDataFrame` with `isinstance(area, gpd.GeoDataFrame)` (`tessellation.py:62`)
- [x] Replace `type(area) == str` with `isinstance(area, str)` (`tessellation.py:64`)
- [x] Replace `type(self.area_gdf) == MultiPolygon` with `isinstance(...)` (`tessellation.py:255, 362`)
- [x] Replace `type(gdf.geometry.iloc[0]) == Polygon` / `MultiPolygon` with `isinstance(...)` (`hexagons.py:27, 41`)
- [x] Replace `type(self.detail_deg) is int` with `isinstance(self.detail_deg, int)` (`roads.py:63`)
- [x] Replace `hasattr(self.poi_dataframe, cat)` with `cat in self.poi_dataframe.columns` (`tessellation.py:81`)
- [x] Replace `hasattr(poi_df, poi_category)` with `poi_category in poi_df.columns` (`poi.py:202`)
- [x] Specify warning category in `warnings.simplefilter("ignore")` — use `FutureWarning` or `ShapelyDeprecationWarning` instead of blanket suppression (`poi.py:68`, `_geo.py:43`, `tessellation.py:415`)
- [x] Remove `E501` from ruff ignore list and fix lines exceeding 88 chars (`pyproject.toml:75`)
- [x] Replace deprecated `affinity="euclidean"` with `metric="euclidean"` in AgglomerativeClustering (`tessellation.py:423`)

---

## Group 3 — Logging

Replace all `print()` statements with Python's `logging` module for proper observability.

**Files:** `tesspy/__init__.py`, `tesspy/tessellation.py`, `tesspy/data/poi.py`, `tesspy/data/roads.py`, `tesspy/_validators.py`

- [x] Add package-level logger in `tesspy/__init__.py`: `logging.getLogger("tesspy")`
- [x] Replace `print()` calls in `tessellation.py` verbose blocks with `logger.info()` (~8 locations)
- [x] Replace `print()` calls in `poi.py` verbose blocks with `logger.info()` / `logger.debug()` (~6 locations)
- [x] Replace `print()` calls in `roads.py` verbose blocks with `logger.info()` (~3 locations)
- [x] Replace `print("MultiPolygon found...")` in `_validators.py:70` with `logger.info()`
- [x] Add `NullHandler` to package logger so library users don't see unwanted output
- [x] Update docstrings to mention logging instead of "print progress information"

---

## Group 4 — Testing Expansion

Add unit tests for untested modules, error-path coverage, and enforce a coverage floor.

**Files:** `tests/unit/` (new files), `tests/conftest.py`, `pyproject.toml`

- [ ] Add `tests/unit/test_hexagons.py` — unit tests for `get_h3_hexagons()` with mock GeoDataFrames (Polygon & MultiPolygon inputs)
- [ ] Add `tests/unit/test_voronoi.py` — unit tests for `voronoi_polygons()` with a small synthetic Voronoi diagram
- [ ] Add `tests/unit/test_clustering.py` — unit tests for `get_hierarchical_clustering_parameter()`
- [ ] Add `tests/unit/test_road_data.py` — unit tests for `RoadData.create_custom_filter()` with various `detail_deg` values
- [ ] Add `tests/unit/test_deprecation_shims.py` — verify that importing from `tesspy.tessellation_functions` and `tesspy.poi_data` emits `DeprecationWarning`
- [ ] Add error-path tests for `POIdata.get_poi_data()` — mock HTTP 429, 504, and non-200 responses
- [ ] Add error-path tests for `Tessellation.__init__()` — invalid area types (int, list, etc.)
- [ ] Add tests for `count_poi_per_tile()` — invalid inputs, empty GeoDataFrame, string city name
- [ ] Add `pytest-timeout` to dev dependencies in `pyproject.toml`
- [ ] Set coverage fail-under threshold (e.g., 70%) in `pyproject.toml` `[tool.coverage.report]`

---

## Group 5 — Documentation

Update outdated docs, fix version strings, and add status badges to README.

**Files:** `docs/source/conf.py`, `docs/source/Contribution.rst`, `README.md`

- [ ] Fix `conf.py` version: replace hardcoded `release = "0.0.1"` with dynamic read from `tesspy/_version.py`
- [ ] Update `conf.py` copyright year from 2022 to 2022–2025
- [ ] Update `Contribution.rst` step 2: replace `python setup.py develop` with `pip install -e ".[dev]"`
- [ ] Update `Contribution.rst` step 7: replace `black` references with `ruff format` and mention pre-commit hooks
- [ ] Update `Contribution.rst` step 2: replace conda-only dependency install with `pip install -e ".[dev]"` option
- [ ] Add CI badge to README (GitHub Actions workflow status)
- [ ] Add PyPI version badge to README
- [ ] Add Python version badge to README
- [ ] Add a "Development" section to README pointing to CONTRIBUTING and this PLAN

---

## Group 6 — CI/CD Hardening

Improve the CI pipeline with security scanning, docs verification, and broader trigger coverage.

**Files:** `.github/workflows/tests_package.yml`, `.github/workflows/docs.yml` (new)

- [ ] Add `develop` branch to push/PR triggers in `tests_package.yml` so CI runs on develop PRs too
- [ ] Add `bandit` security scanning step (or `ruff` S rules) to lint job
- [ ] Add `pip-audit` dependency vulnerability scanning step
- [ ] Add Sphinx docs build verification workflow (`.github/workflows/docs.yml`) — build docs and fail on warnings
- [ ] Add coverage `--fail-under` flag to pytest step in CI
- [ ] Add `pytest-timeout` to the test runner with a global timeout (e.g., 300s for unit tests)
- [ ] Pin action versions with full SHA hashes instead of tags for supply-chain security

---

## Group 7 — API Resilience & Robustness

Add retry logic, timeout handling, and input validation to protect against transient API failures.

**Files:** `tesspy/data/poi.py`, `tesspy/data/roads.py`, `tesspy/data/_geo.py`, `tesspy/_validators.py`

- [ ] Add configurable `requests.Session` with retry adapter (urllib3 `Retry`) for Overpass API calls in `poi.py`
- [ ] Add `requests` timeout parameter to `get_poi_data()` HTTP call (currently no socket-level timeout)
- [ ] Add retry/backoff logic for osmnx calls in `roads.py` (network failures, rate limiting)
- [ ] Add geometry validity check in `_validators.py` — detect self-intersecting polygons with `shapely.validation.make_valid()`
- [ ] Add explicit check for empty/None geometries before spatial operations
- [ ] Improve generic "Bad Request!" error message in `poi.py:131` — include status code and response body excerpt
- [ ] Add `ConnectionError` / `Timeout` exception handling around network calls with helpful user-facing messages

---

## Group 8 — Release Preparation (v0.2.0)

Final steps to ship the v0.2.0 release.

**Files:** `tesspy/_version.py`, `CHANGELOG.md`, `tesspy/tessellation_functions.py`, `tesspy/poi_data.py`

- [ ] Bump version in `_version.py` from `"0.1.2"` to `"0.2.0"`
- [ ] Finalize CHANGELOG.md — move all items from `[Unreleased]` to `[0.2.0] — YYYY-MM-DD`
- [ ] Add deprecation timeline to shim module docstrings: "Will be removed in v0.3.0"
- [ ] Run full test suite (unit + integration) and confirm all pass
- [ ] Verify `python -m build` produces correct wheel and sdist
- [ ] Verify OIDC trusted publishing is configured in PyPI project settings
- [ ] Tag release `v0.2.0` and create GitHub Release with notes from CHANGELOG
- [ ] Publish to PyPI
- [ ] Update conda-forge feedstock recipe (if applicable)

---

## Suggested execution order

1. **Group 2** (Code Quality) — quick wins, no new files
2. **Group 3** (Logging) — builds on Group 2's cleaner code
3. **Group 5** (Documentation) — can be done in parallel with Groups 2–3
4. **Group 4** (Testing) — write tests against the cleaned-up code
5. **Group 6** (CI/CD) — enforce the new tests and quality gates
6. **Group 7** (API Resilience) — deeper changes, needs good test coverage first
7. **Group 8** (Release) — final step once everything above is done
