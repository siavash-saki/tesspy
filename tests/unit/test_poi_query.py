"""
Unit tests for POI category validation and road data filter generation
(no network required).
"""

import geopandas as gpd
import pytest
from shapely.geometry import Point

from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData


def test_invalid_category_raises_on_get_poi_data(soho_polygon_gdf, monkeypatch):
    """An invalid POI category should raise ValueError in get_poi_data."""
    poi_data = POIdata(
        soho_polygon_gdf,
        poi_categories=["not_a_real_category"],
        timeout=60,
        verbose=False,
    )
    with pytest.raises(ValueError, match="not a valid POI primary category"):
        poi_data.get_poi_data()


def test_valid_category_calls_osmnx(soho_polygon_gdf, monkeypatch):
    """Valid categories should pass validation and call osmnx."""
    import tesspy.data.poi as poi_module

    fake_gdf = gpd.GeoDataFrame(
        {"amenity": ["cafe"]},
        geometry=[Point(-0.13, 51.51)],
        crs="EPSG:4326",
    )
    monkeypatch.setattr(
        poi_module.ox, "features_from_polygon", lambda *a, **kw: fake_gdf
    )

    poi_data = POIdata(
        soho_polygon_gdf,
        poi_categories=["amenity"],
        timeout=60,
        verbose=False,
    )
    result = poi_data.get_poi_data()
    assert "center_longitude" in result.columns
    assert "center_latitude" in result.columns
    assert "amenity" in result.columns
    assert len(result) == 1


def test_create_custom_filter_all_types(soho_polygon_gdf):
    road_data = RoadData(soho_polygon_gdf, detail_deg=None)
    cf = road_data.create_custom_filter()

    assert isinstance(cf, str)
    assert "motorway" in cf
    assert "trunk" in cf
    assert "residential" in cf
    assert "|" in cf
    assert "amenity" not in cf
    assert "public_transport" not in cf


def test_create_custom_filter_limited_types(soho_polygon_gdf):
    road_data = RoadData(soho_polygon_gdf, detail_deg=3)
    cf = road_data.create_custom_filter()

    assert "motorway" in cf
    assert "trunk" in cf
    assert "primary" in cf
    assert "secondary" not in cf


def test_create_custom_filter_invalid_detail_deg(soho_polygon_gdf):
    road_data = RoadData(soho_polygon_gdf, detail_deg="invalid")
    with pytest.raises(ValueError, match="detail_deg must be None or an int"):
        road_data.create_custom_filter()
