"""
Unit tests for POIdata.get_poi_data error paths (no network required).
Uses monkeypatch to stub HTTP responses.
"""

import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from tesspy.data.poi import POIdata


@pytest.fixture()
def area_gdf():
    """A small polygon for constructing POIdata."""
    poly = Polygon([
        (13.38, 52.51),
        (13.42, 52.51),
        (13.42, 52.53),
        (13.38, 52.53),
        (13.38, 52.51),
    ])
    return gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")


class _FakeResponse:
    """Minimal mock for requests.Response."""

    def __init__(self, status_code, text=""):
        self.status_code = status_code
        self.text = text


def test_http_429_raises_runtime_error(area_gdf, monkeypatch):
    """HTTP 429 should raise RuntimeError with '429' in the message."""
    monkeypatch.setattr(
        "tesspy.data.poi.requests.get",
        lambda **kwargs: _FakeResponse(429),
    )
    poi = POIdata(area_gdf, ["amenity"], timeout=10, verbose=False)
    with pytest.raises(RuntimeError, match="429"):
        poi.get_poi_data()


def test_http_504_raises_runtime_error(area_gdf, monkeypatch):
    """HTTP 504 should raise RuntimeError with '504' in the message."""
    monkeypatch.setattr(
        "tesspy.data.poi.requests.get",
        lambda **kwargs: _FakeResponse(504),
    )
    poi = POIdata(area_gdf, ["amenity"], timeout=10, verbose=False)
    with pytest.raises(RuntimeError, match="504"):
        poi.get_poi_data()


def test_http_other_error_raises_runtime_error(area_gdf, monkeypatch):
    """Non-200/429/504 status should raise RuntimeError with 'Bad Request'."""
    monkeypatch.setattr(
        "tesspy.data.poi.requests.get",
        lambda **kwargs: _FakeResponse(400, text="error details here"),
    )
    poi = POIdata(area_gdf, ["amenity"], timeout=10, verbose=False)
    with pytest.raises(RuntimeError, match="Bad Request"):
        poi.get_poi_data()


def test_invalid_category_raises_value_error(area_gdf):
    """An invalid POI category should raise ValueError."""
    poi = POIdata(area_gdf, ["not_a_real_category"], timeout=10, verbose=False)
    with pytest.raises(ValueError, match="not a valid POI primary category"):
        poi.get_poi_data()


def test_empty_response_raises_value_error(area_gdf, monkeypatch):
    """HTTP 200 with no elements should raise ValueError."""
    import json

    empty_resp = json.dumps({"elements": []})
    monkeypatch.setattr(
        "tesspy.data.poi.requests.get",
        lambda **kwargs: _FakeResponse(200, text=empty_resp),
    )
    poi = POIdata(area_gdf, ["amenity"], timeout=10, verbose=False)
    with pytest.raises(ValueError, match="No POI data found"):
        poi.get_poi_data()
