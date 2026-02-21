"""
Unit tests for Tessellation.__init__ error paths (no network required).
"""

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, Polygon

from tesspy.tessellation import Tessellation


@pytest.fixture()
def valid_gdf():
    """A valid single-Polygon GeoDataFrame with CRS."""
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    return gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")


def test_valid_geodataframe_accepted(valid_gdf):
    """Tessellation should accept a valid GeoDataFrame."""
    t = Tessellation(valid_gdf)
    assert t.area_gdf is not None
    assert isinstance(t.poi_dataframe, pd.DataFrame)
    assert t.available_poi_categories == []
    assert t.queried_highway_types == []


def test_invalid_type_int_raises():
    """Passing an int should raise TypeError."""
    with pytest.raises(TypeError, match="area must be a GeoDataFrame or a string"):
        Tessellation(42)


def test_invalid_type_list_raises():
    """Passing a list should raise TypeError."""
    with pytest.raises(TypeError, match="area must be a GeoDataFrame or a string"):
        Tessellation([1, 2, 3])


def test_invalid_type_none_raises():
    """Passing None should raise TypeError."""
    with pytest.raises(TypeError, match="area must be a GeoDataFrame or a string"):
        Tessellation(None)


def test_empty_geodataframe_raises():
    """A GeoDataFrame with no rows should raise ValueError."""
    empty_gdf = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    with pytest.raises(ValueError, match="exactly one geometry element"):
        Tessellation(empty_gdf)


def test_multiple_rows_raises():
    """A GeoDataFrame with >1 rows should raise ValueError."""
    poly1 = Polygon([(0, 0), (1, 0), (1, 1)])
    poly2 = Polygon([(2, 2), (3, 2), (3, 3)])
    gdf = gpd.GeoDataFrame(geometry=[poly1, poly2], crs="EPSG:4326")
    with pytest.raises(ValueError, match="exactly one geometry element"):
        Tessellation(gdf)


def test_wrong_geometry_type_raises():
    """A GeoDataFrame with a Point geometry should raise TypeError."""
    gdf = gpd.GeoDataFrame(geometry=[Point(0, 0)], crs="EPSG:4326")
    with pytest.raises(TypeError, match="Polygon or MultiPolygon"):
        Tessellation(gdf)


def test_line_geometry_raises():
    """A GeoDataFrame with a LineString should raise TypeError."""
    gdf = gpd.GeoDataFrame(geometry=[LineString([(0, 0), (1, 1)])], crs="EPSG:4326")
    with pytest.raises(TypeError, match="Polygon or MultiPolygon"):
        Tessellation(gdf)


def test_no_crs_raises():
    """A GeoDataFrame without CRS should raise ValueError."""
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    gdf = gpd.GeoDataFrame(geometry=[poly])
    with pytest.raises(ValueError, match="must have a CRS"):
        Tessellation(gdf)
