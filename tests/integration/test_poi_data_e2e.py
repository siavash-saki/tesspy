"""
Integration tests for POI and road data retrieval (requires OSM API access).
"""

import geopandas as gpd
import pytest

from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData
from tesspy.tessellation import Tessellation

from tests.conftest import call_with_osm_retry


@pytest.mark.integration
@pytest.mark.slow
def test_get_poi_data_city_1():
    city1 = Tessellation("Innenstadt, Frankfurt am Main")
    poi_data = POIdata(
        city1.get_polygon(),
        poi_categories=["public_transport"],
        timeout=60,
        verbose=False,
    )

    ffm_data = call_with_osm_retry(poi_data.get_poi_data)

    assert len(ffm_data) > 0
    assert hasattr(ffm_data, "geometry")
    assert hasattr(ffm_data, "public_transport")
    assert hasattr(ffm_data, "center_latitude")
    assert hasattr(ffm_data, "center_longitude")
    assert hasattr(ffm_data, "tags")


@pytest.mark.integration
@pytest.mark.slow
def test_get_poi_data_city_2():
    city2 = Tessellation("Downtown San Diego")
    poi_data_2 = POIdata(
        city2.get_polygon(),
        poi_categories=["public_transport"],
        timeout=60,
        verbose=False,
    )

    dsd_data = call_with_osm_retry(poi_data_2.get_poi_data)

    assert len(dsd_data) > 0
    assert hasattr(dsd_data, "geometry")
    assert hasattr(dsd_data, "public_transport")
    assert hasattr(dsd_data, "center_latitude")
    assert hasattr(dsd_data, "center_longitude")
    assert hasattr(dsd_data, "tags")


@pytest.mark.integration
@pytest.mark.slow
def test_get_road_network_data():
    city = Tessellation("SOHO, London").get_polygon()
    road_data = RoadData(city).get_road_network()

    assert type(road_data) == gpd.GeoDataFrame
    assert len(road_data) > 100
    assert hasattr(road_data, "geometry")
    assert hasattr(road_data, "highway")
    assert hasattr(road_data, "osmid")


@pytest.mark.integration
@pytest.mark.slow
def test_get_road_network_split():
    city = Tessellation("Innenstadt, Frankfurt am Main").get_polygon()
    split_road_data = RoadData(city, split_roads=True, verbose=True).get_road_network()
    road_data = RoadData(city, split_roads=False, verbose=True).get_road_network()

    assert len(split_road_data) > 0
    assert len(road_data) > 0
    assert len(split_road_data) >= len(road_data)
