from tesspy.poi_data import *
from tesspy.tessellation import *
from tesspy.tessellation_functions import *


def test_count_poi():
    city = Tessellation("Berlin")
    city_polygon = city.get_polygon()
    city_poi = POIdata(city_polygon, poi_categories=["leisure"], timeout=60, verbose=False).get_poi_data()
    points_geom = city_poi[['center_longitude', 'center_latitude']]\
        .apply(lambda p: Point(p['center_longitude'], p['center_latitude']), axis=1)
    tess_data = gpd.GeoDataFrame(geometry=points_geom,
                                 crs='EPSG:4326')
    city_squares = city.squares(14)
    count1 = count_poi(city_squares, tess_data)
    assert hasattr(count1, "count")
    assert len(count1) > 0


def test_get_adaptive_squares():
    city = Tessellation("San Diego")
    city_poi = POIdata(city.get_polygon(), poi_categories=["leisure"], timeout=60, verbose=False).get_poi_data()
    points_geom = city_poi[['center_longitude', 'center_latitude']] \
        .apply(lambda p: Point(p['center_longitude'], p['center_latitude']), axis=1)
    tess_data = gpd.GeoDataFrame(geometry=points_geom,
                                 crs='EPSG:4326')
    city_squares = city.squares(14)
    count1 = count_poi(city_squares, tess_data)

    assert len(get_adaptive_squares(count1, 100)) > len(city_squares)
    assert len(get_adaptive_squares(count1, len(tess_data))) == len(city_squares)


def test_split_linestring():
    city = Tessellation("Lille")
    road_data = RoadData(city.get_polygon(), split_roads=False, verbose=True).get_road_network()
    road_data2 = RoadData(city.get_polygon(), split_roads=True, verbose=False).get_road_network()
    split_road_data = split_linestring(road_data)
    assert len(split_road_data) > len(road_data)
    assert len(split_road_data) == road_data2
