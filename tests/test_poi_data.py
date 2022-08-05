from tesspy.tessellation_functions import *
from tesspy.poi_data import *
from tesspy.tessellation import *
import geopandas as gpd
from time import sleep


def test_osm_primary_features():
    primary_categories = [
        "aerialway",
        "aeroway",
        "amenity",
        "barrier",
        "boundary",
        "building",
        "craft",
        "emergency",
        "geological",
        "healthcare",
        "highway",
        "historic",
        "landuse",
        "leisure",
        "man_made",
        "military",
        "natural",
        "office",
        "place",
        "power",
        "public_transport",
        "railway",
        "route",
        "shop",
        "sport",
        "telecom",
        "tourism",
        "water",
        "waterway",
    ]
    city = Tessellation("Berlin")
    poi_data = POIdata(
        city.get_polygon(), poi_categories=["landuse"], timeout=60, verbose=False
    )
    ls_osm_features = poi_data.osm_primary_features()
    assert type(ls_osm_features) == list
    assert len(ls_osm_features) > 0
    assert primary_categories == ls_osm_features
    assert "house" not in ls_osm_features


def test_create_overpass_query_string():
    city = Tessellation("Frankfurt am Main")
    poi_data = POIdata(
        city.get_polygon(),
        poi_categories=["public_transport"],
        timeout=60,
        verbose=False,
    )
    poi_data2 = POIdata(
        city.get_polygon(),
        poi_categories=["public_transport", "tourism"],
        timeout=60,
        verbose=False,
    )
    poi_data3 = POIdata(
        city.get_polygon(),
        poi_categories=["public_transport", "shop", "waterway"],
        timeout=60,
        verbose=False,
    )
    query_str = poi_data.create_overpass_query_string()
    assert type(query_str) == str
    assert len(query_str) > 0
    assert "timeout:60" in query_str
    assert "out:json" in query_str
    assert "node[public_transport]" in query_str
    assert "node[tourism]" not in query_str
    assert (
        "node[public_transport]" in poi_data2.create_overpass_query_string()
        and "node[tourism]" in poi_data2.create_overpass_query_string()
    )
    assert "way[public_transport]" in poi_data3.create_overpass_query_string()


def test_get_poi_data_city_1():
    # City 1
    city1 = Tessellation("Innenstadt, Frankfurt am Main")
    poi_data = POIdata(
        city1.get_polygon(),
        poi_categories=["public_transport"],
        timeout=60,
        verbose=False,
    )

    while True:
        try:
            ffm_data = poi_data.get_poi_data()
            break
        except RuntimeError:
            sleep(5 * 60)
            continue

    assert len(ffm_data) > 0
    assert hasattr(ffm_data, "geometry")
    assert hasattr(ffm_data, "public_transport")
    assert hasattr(ffm_data, "center_latitude")
    assert hasattr(ffm_data, "center_longitude")
    assert hasattr(ffm_data, "tags")


def test_get_poi_data_city_2():
    sleep(10)
    # City 2
    city2 = Tessellation("Downtown San Diego")
    poi_data_2 = POIdata(
        city2.get_polygon(),
        poi_categories=["public_transport"],
        timeout=60,
        verbose=False,
    )

    while True:
        try:
            dsd_data = poi_data_2.get_poi_data()
            break
        except RuntimeError:
            sleep(5 * 60)
            continue

    assert len(dsd_data) > 0
    assert hasattr(dsd_data, "geometry")
    assert hasattr(dsd_data, "public_transport")
    assert hasattr(dsd_data, "center_latitude")
    assert hasattr(dsd_data, "center_longitude")
    assert hasattr(dsd_data, "tags")


def test_osm_highway_types():
    city = Tessellation("Key West")
    road_data = RoadData(city)
    ls_highway_types = road_data.osm_highway_types()

    assert type(ls_highway_types) == list
    assert len(ls_highway_types) > 0
    assert "primary" in ls_highway_types
    assert "secondary" in ls_highway_types
    assert "residential" in ls_highway_types
    assert "street" not in ls_highway_types


def test_create_custom_filter():
    city = Tessellation("Johannesburg")
    cf = RoadData(city).create_custom_filter()
    assert type(cf) == str
    assert len(cf) > 0
    assert "motorway" in cf
    assert "trunk" in cf
    assert "residential" in cf
    assert "|" in cf
    assert "amenity" not in cf
    assert "public_transport" not in cf


def test_get_road_network_data():
    city = Tessellation("SOHO, London").get_polygon()
    sleep(60)
    road_data = RoadData(city).get_road_network()
    assert type(road_data) == gpd.GeoDataFrame
    assert len(road_data) > 100
    assert hasattr(road_data, "geometry")
    assert hasattr(road_data, "highway")
    assert hasattr(road_data, "osmid")


def test_get_road_network_split():
    city = Tessellation("Innenstadt, Frankfurt am Main").get_polygon()
    sleep(60)
    split_road_data = RoadData(city, split_roads=True, verbose=True).get_road_network()
    sleep(60)
    road_data = RoadData(city, split_roads=False, verbose=True).get_road_network()
    assert len(split_road_data) > 0
    assert len(road_data) > 0
    assert len(split_road_data) >= len(road_data)
