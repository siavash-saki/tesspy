"""
Unit tests for OSM constant lists (no network required).
"""

from tesspy._constants import OSM_HIGHWAY_TYPES, OSM_PRIMARY_FEATURES
from tesspy.tessellation import Tessellation
from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData


def test_osm_primary_features_content():
    expected = [
        "aerialway", "aeroway", "amenity", "barrier", "boundary", "building",
        "craft", "emergency", "geological", "healthcare", "highway", "historic",
        "landuse", "leisure", "man_made", "military", "natural", "office",
        "place", "power", "public_transport", "railway", "route", "shop",
        "sport", "telecom", "tourism", "water", "waterway",
    ]
    assert OSM_PRIMARY_FEATURES == expected
    assert "house" not in OSM_PRIMARY_FEATURES


def test_osm_primary_features_type():
    assert isinstance(OSM_PRIMARY_FEATURES, list)
    assert len(OSM_PRIMARY_FEATURES) == 29
    assert all(isinstance(f, str) for f in OSM_PRIMARY_FEATURES)


def test_osm_highway_types_content():
    assert "motorway" in OSM_HIGHWAY_TYPES
    assert "primary" in OSM_HIGHWAY_TYPES
    assert "secondary" in OSM_HIGHWAY_TYPES
    assert "residential" in OSM_HIGHWAY_TYPES
    assert "street" not in OSM_HIGHWAY_TYPES


def test_osm_highway_types_type():
    assert isinstance(OSM_HIGHWAY_TYPES, list)
    assert len(OSM_HIGHWAY_TYPES) == 19
    assert all(isinstance(t, str) for t in OSM_HIGHWAY_TYPES)


def test_tessellation_static_methods_delegate_to_constants():
    """Verify the Tessellation static methods return the shared constants."""
    assert Tessellation.osm_primary_features() is OSM_PRIMARY_FEATURES
    assert Tessellation.osm_highway_types() is OSM_HIGHWAY_TYPES


def test_poi_data_static_method_delegates():
    assert POIdata.osm_primary_features() is OSM_PRIMARY_FEATURES


def test_road_data_static_method_delegates():
    assert RoadData.osm_highway_types() is OSM_HIGHWAY_TYPES
