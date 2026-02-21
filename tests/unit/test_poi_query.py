"""
Unit tests for Overpass API query string generation (no network required).
Uses hardcoded polygon fixtures instead of geocoding calls.
"""

import pytest

from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData


def test_create_overpass_query_string_single_category(soho_polygon_gdf):
    poi_data = POIdata(
        soho_polygon_gdf,
        poi_categories=["public_transport"],
        timeout=60,
        verbose=False,
    )
    query = poi_data.create_overpass_query_string()

    assert isinstance(query, str)
    assert len(query) > 0
    assert "timeout:60" in query
    assert "out:json" in query
    assert "node[public_transport]" in query
    assert "node[tourism]" not in query


def test_create_overpass_query_string_multiple_categories(soho_polygon_gdf):
    poi_data = POIdata(
        soho_polygon_gdf,
        poi_categories=["public_transport", "tourism"],
        timeout=60,
        verbose=False,
    )
    query = poi_data.create_overpass_query_string()

    assert "node[public_transport]" in query
    assert "node[tourism]" in query


def test_create_overpass_query_string_three_categories(soho_polygon_gdf):
    poi_data = POIdata(
        soho_polygon_gdf,
        poi_categories=["public_transport", "shop", "waterway"],
        timeout=60,
        verbose=False,
    )
    query = poi_data.create_overpass_query_string()

    assert "way[public_transport]" in query
    assert "way[shop]" in query
    assert "way[waterway]" in query


def test_create_overpass_query_invalid_category(soho_polygon_gdf):
    poi_data = POIdata(
        soho_polygon_gdf,
        poi_categories=["not_a_real_category"],
        timeout=60,
        verbose=False,
    )
    with pytest.raises(ValueError, match="not a valid POI primary category"):
        poi_data.create_overpass_query_string()


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
