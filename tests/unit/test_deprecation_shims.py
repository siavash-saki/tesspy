"""
Unit tests for deprecation shim modules (no network required).
Verify that importing from legacy module paths emits DeprecationWarning.
"""

import importlib
import sys
import warnings

import pytest


def _fresh_import(module_name: str):
    """Force a fresh import of a module by removing it from sys.modules."""
    sys.modules.pop(module_name, None)
    return importlib.import_module(module_name)


def test_tessellation_functions_emits_deprecation_warning():
    """Importing tesspy.tessellation_functions should emit DeprecationWarning."""
    sys.modules.pop("tesspy.tessellation_functions", None)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _fresh_import("tesspy.tessellation_functions")

    deprecation_warnings = [
        w for w in caught if issubclass(w.category, DeprecationWarning)
    ]
    assert len(deprecation_warnings) >= 1
    assert "tessellation_functions" in str(deprecation_warnings[0].message)


def test_poi_data_emits_deprecation_warning():
    """Importing tesspy.poi_data should emit DeprecationWarning."""
    sys.modules.pop("tesspy.poi_data", None)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _fresh_import("tesspy.poi_data")

    deprecation_warnings = [
        w for w in caught if issubclass(w.category, DeprecationWarning)
    ]
    assert len(deprecation_warnings) >= 1
    assert "poi_data" in str(deprecation_warnings[0].message)


def test_tessellation_functions_re_exports_names():
    """The shim should re-export all expected names."""
    sys.modules.pop("tesspy.tessellation_functions", None)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        mod = _fresh_import("tesspy.tessellation_functions")

    expected = [
        "count_poi",
        "get_squares_polyfill",
        "get_adaptive_squares",
        "get_h3_hexagons",
        "voronoi_polygons",
        "create_city_blocks",
        "merge_city_blocks",
        "split_linestring",
        "explode",
        "create_blocks",
        "get_rest_polygon",
        "get_hierarchical_clustering_parameter",
    ]
    for name in expected:
        assert hasattr(mod, name), f"Missing re-export: {name}"


def test_poi_data_re_exports_names():
    """The shim should re-export all expected names."""
    sys.modules.pop("tesspy.poi_data", None)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        mod = _fresh_import("tesspy.poi_data")

    for name in ["POIdata", "RoadData", "geom_ceil", "geom_floor"]:
        assert hasattr(mod, name), f"Missing re-export: {name}"
