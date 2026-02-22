"""
Unit tests for POIdata.get_poi_data error paths (no network required).
Uses monkeypatch to stub osmnx calls.
"""

import geopandas as gpd
import pytest
from shapely.geometry import Point, Polygon

import tesspy.data.poi as poi_module
from tesspy.data.poi import POIdata


@pytest.fixture()
def area_gdf():
    """A small polygon for constructing POIdata."""
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


def test_invalid_category_raises_value_error(area_gdf):
    """An invalid POI category should raise ValueError."""
    poi = POIdata(area_gdf, ["not_a_real_category"], timeout=10, verbose=False)
    with pytest.raises(ValueError, match="not a valid POI primary category"):
        poi.get_poi_data()


def test_empty_response_raises_value_error(area_gdf, monkeypatch):
    """osmnx returning an empty GeoDataFrame should raise ValueError."""
    empty_gdf = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    monkeypatch.setattr(
        poi_module.ox,
        "features_from_polygon",
        lambda *args, **kwargs: empty_gdf,
    )
    poi = POIdata(area_gdf, ["amenity"], timeout=10, verbose=False)
    with pytest.raises(ValueError, match="No POI data found"):
        poi.get_poi_data()


def test_osmnx_exception_propagates(area_gdf, monkeypatch):
    """Network errors from osmnx should propagate to the caller."""

    def raise_error(*args, **kwargs):
        raise RuntimeError("Overpass server error")

    monkeypatch.setattr(
        poi_module.ox,
        "features_from_polygon",
        raise_error,
    )
    poi = POIdata(area_gdf, ["amenity"], timeout=10, verbose=False)
    with pytest.raises(RuntimeError, match="Overpass server error"):
        poi.get_poi_data()


def test_missing_category_column_returns_false(area_gdf, monkeypatch):
    """If a queried category has no results, the column should be all False."""
    fake_gdf = gpd.GeoDataFrame(
        {"amenity": ["cafe"]},
        geometry=[Point(13.40, 52.52)],
        crs="EPSG:4326",
    )
    monkeypatch.setattr(
        poi_module.ox,
        "features_from_polygon",
        lambda *args, **kwargs: fake_gdf,
    )
    poi = POIdata(area_gdf, ["amenity", "building"], timeout=10, verbose=False)
    result = poi.get_poi_data()

    assert result["amenity"].iloc[0]
    assert not result["building"].iloc[0]
