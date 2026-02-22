"""
Unit tests for hexagon tessellation functions (no network required).
"""

import geopandas as gpd
import pytest
from shapely.geometry import MultiPolygon, Polygon

from tesspy.methods.hexagons import get_h3_hexagons


@pytest.fixture()
def small_polygon_gdf():
    """A small polygon in central Berlin large enough for H3 resolution 7."""
    poly = Polygon(
        [
            (13.38, 52.51),
            (13.42, 52.51),
            (13.42, 52.53),
            (13.38, 52.53),
            (13.38, 52.51),
        ]
    )
    return gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")


@pytest.fixture()
def multi_polygon_gdf():
    """A MultiPolygon consisting of two separate squares, large enough for H3."""
    poly1 = Polygon(
        [
            (13.30, 52.48),
            (13.40, 52.48),
            (13.40, 52.54),
            (13.30, 52.54),
            (13.30, 52.48),
        ]
    )
    poly2 = Polygon(
        [
            (13.45, 52.48),
            (13.55, 52.48),
            (13.55, 52.54),
            (13.45, 52.54),
            (13.45, 52.48),
        ]
    )
    mp = MultiPolygon([poly1, poly2])
    return gpd.GeoDataFrame(geometry=[mp], crs="EPSG:4326")


def test_polygon_returns_geodataframe(small_polygon_gdf):
    """get_h3_hexagons returns a GeoDataFrame with hexagons for a Polygon input."""
    result = get_h3_hexagons(small_polygon_gdf, resolution=7)

    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) > 0
    assert "geometry" in result.columns
    assert result.crs is not None


def test_polygon_geometries_are_polygons(small_polygon_gdf):
    """All returned geometries should be Polygon instances."""
    result = get_h3_hexagons(small_polygon_gdf, resolution=7)

    for geom in result.geometry:
        assert isinstance(geom, Polygon)


def test_higher_resolution_more_hexagons(small_polygon_gdf):
    """Higher H3 resolution should produce more (smaller) hexagons."""
    low_res = get_h3_hexagons(small_polygon_gdf, resolution=5)
    high_res = get_h3_hexagons(small_polygon_gdf, resolution=7)

    assert len(high_res) > len(low_res)


def test_multipolygon_returns_geodataframe(multi_polygon_gdf):
    """get_h3_hexagons handles MultiPolygon input."""
    result = get_h3_hexagons(multi_polygon_gdf, resolution=7)

    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) > 0
    assert "geometry" in result.columns


def test_index_contains_h3_ids(small_polygon_gdf):
    """The returned GeoDataFrame index should contain H3 hex ID strings."""
    result = get_h3_hexagons(small_polygon_gdf, resolution=7)

    for idx in result.index:
        assert isinstance(idx, str)
        assert len(idx) > 0
