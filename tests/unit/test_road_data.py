"""
Unit tests for RoadData.create_custom_filter (no network required).
"""

import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from tesspy._constants import OSM_HIGHWAY_TYPES
from tesspy.data.roads import RoadData


@pytest.fixture()
def dummy_area():
    """Minimal GeoDataFrame for constructing RoadData."""
    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    return gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")


def test_none_detail_deg_uses_all_types(dummy_area):
    """detail_deg=None should include all 19 highway types in the filter."""
    rd = RoadData(dummy_area, detail_deg=None)
    filt = rd.create_custom_filter()

    for hw_type in OSM_HIGHWAY_TYPES:
        assert hw_type in filt

    assert filt.startswith("['highway'~'")
    assert filt.endswith("']")


def test_detail_deg_slices_types(dummy_area):
    """detail_deg=3 should include only the first 3 highway types."""
    rd = RoadData(dummy_area, detail_deg=3)
    filt = rd.create_custom_filter()

    for hw_type in OSM_HIGHWAY_TYPES[:3]:
        assert hw_type in filt

    # Types beyond the slice should not appear
    for hw_type in OSM_HIGHWAY_TYPES[3:]:
        assert hw_type not in filt


def test_detail_deg_invalid_type_raises(dummy_area):
    """Non-int, non-None detail_deg should raise ValueError."""
    rd = RoadData(dummy_area, detail_deg="high")
    with pytest.raises(ValueError, match="detail_deg must be None or an int"):
        rd.create_custom_filter()


def test_filter_format(dummy_area):
    """Filter should follow osmnx format: ['highway'~'type1|type2|...']."""
    rd = RoadData(dummy_area, detail_deg=2)
    filt = rd.create_custom_filter()

    assert filt == "['highway'~'motorway|trunk']"
