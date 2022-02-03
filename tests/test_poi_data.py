from tesspy.tessellation_functions import *
from tesspy.poi_data import *
from tesspy.tessellation import *
from time import sleep


def test_osm_primary_features():
    primary_categories = ['aerialway',
                          'aeroway',
                          'amenity',
                          'barrier',
                          'boundary',
                          'building',
                          'craft',
                          'emergency',
                          'geological',
                          'healthcare',
                          'highway',
                          'historic',
                          'landuse',
                          'leisure',
                          'man_made',
                          'military',
                          'natural',
                          'office',
                          'place',
                          'power',
                          'public_transport',
                          'railway',
                          'route',
                          'shop',
                          'sport',
                          'telecom',
                          'tourism',
                          'water',
                          'waterway']
    city = Tessellation("Frankfurt am Main")
    poi_data = POIdata(city.get_polygon(), poi_categories=["amenity"], timeout=60, verbose=False)
    assert primary_categories == poi_data.osm_primary_features()
    assert "house" not in poi_data.osm_primary_features()


def test_create_overpass_query_string():
    city = Tessellation("Frankfurt am Main")
    poi_data = POIdata(city.get_polygon(), poi_categories=["amenity"], timeout=60, verbose=False)
    poi_data2 = POIdata(city.get_polygon(), poi_categories=["amenity", "tourism"], timeout=60, verbose=False)
    poi_data3 = POIdata(city.get_polygon(), poi_categories=["public_transport", "shop", "waterway"], timeout=60,
                        verbose=False)
    assert "node[amenity]" in poi_data.create_overpass_query_string()
    assert "node[amenity]" in poi_data2.create_overpass_query_string() \
           and "node[tourism]" in poi_data2.create_overpass_query_string()
    assert "way[public_transport]" in poi_data3.create_overpass_query_string()


def test_get_poi_data():
    city1 = Tessellation("Frankfurt am Main")
    poi_data = POIdata(city1.get_polygon(), poi_categories=["amenity"], timeout=60, verbose=False)
    assert len(poi_data.get_poi_data()) > 0

    sleep(60)

    city2 = Tessellation("San Diego")
    poi_data2 = POIdata(city2.get_polygon(), poi_categories=["public_transport"], timeout=60, verbose=False)
    print(poi_data2.get_poi_data().columns)
    assert hasattr(poi_data2.get_poi_data(), "geometry")


def test_osm_highway_types():
    city = Tessellation("Key West")
    road_data = RoadData(city)
    assert "primary" in road_data.osm_highway_types()
    assert "street" not in road_data.osm_highway_types()


def test_create_custom_filter():
    city = Tessellation("Johannesburg")
    cf = RoadData(city).create_custom_filter()
    assert "motorway" in cf
    assert "|" in cf
    assert "amenity" not in cf
    assert "residential" in cf


def test_get_road_network():
    city = Tessellation("Johannesburg").get_polygon()
    split_road_data = RoadData(city, split_roads=True, verbose=True).get_road_network()
    road_data = RoadData(city, split_roads=False, verbose=True).get_road_network()
    assert len(split_road_data) > 0
    assert len(road_data) > 0
    assert len(split_road_data) > len(road_data)

