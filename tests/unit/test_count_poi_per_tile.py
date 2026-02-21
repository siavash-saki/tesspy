"""
Unit tests for count_poi_per_tile input validation (no network required).
"""

import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from tesspy.data._geo import count_poi_per_tile


@pytest.fixture()
def simple_gdf():
    """A minimal tessellation GeoDataFrame with one tile."""
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    return gpd.GeoDataFrame(
        {"tile_id": ["tile0"], "geometry": [poly]}, crs="EPSG:4326"
    )


@pytest.fixture()
def area_gdf():
    """A valid area GeoDataFrame."""
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    return gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")


def test_invalid_city_type_raises(simple_gdf):
    """Passing an int as city should raise TypeError."""
    with pytest.raises(TypeError, match="city must be a GeoDataFrame or a city name"):
        count_poi_per_tile(city=123, gdf=simple_gdf)


def test_empty_gdf_raises(area_gdf):
    """An empty tessellation GeoDataFrame should raise ValueError."""
    empty_gdf = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    with pytest.raises(ValueError, match="at least one tile"):
        count_poi_per_tile(city=area_gdf, gdf=empty_gdf)


def test_invalid_poi_categories_type_raises(area_gdf, simple_gdf):
    """Non-string, non-list poi_categories should raise ValueError."""
    with pytest.raises(ValueError, match="poi_categories must be a string or list"):
        count_poi_per_tile(
            city=area_gdf, gdf=simple_gdf, poi_categories=42
        )


def test_string_poi_category_accepted(area_gdf, simple_gdf, monkeypatch):
    """A string poi_category should be wrapped in a list and passed through."""
    # We only want to test the validation logic, not the actual API call,
    # so we stub POIdata.get_poi_data to raise immediately after validation passes.
    def fake_get_poi_data(self):
        raise ValueError("stubbed — validation passed")

    monkeypatch.setattr(
        "tesspy.data._geo.POIdata.get_poi_data", fake_get_poi_data
    )
    # If validation of poi_categories fails, we'd get a different error.
    # We expect the stub's ValueError, confirming the string was accepted.
    with pytest.raises(ValueError, match="stubbed"):
        count_poi_per_tile(
            city=area_gdf, gdf=simple_gdf, poi_categories="amenity"
        )
