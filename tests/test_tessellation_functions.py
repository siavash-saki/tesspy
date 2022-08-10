from tesspy.poi_data import *
from tesspy.tessellation import *
from tesspy.tessellation_functions import *
import pandas as pd
import geopandas as gpd

tess_data = gpd.read_file("tests/leisure_poi_Mitte_Berlin.geojson")
city = Tessellation("Mitte, Berlin")
city_polygon = city.get_polygon()


def test_count_poi():
    city_squares = get_squares_polyfill(city_polygon, 14)
    count1 = count_poi(city_squares, tess_data)

    assert hasattr(count1, "count")
    assert hasattr(count1, "children_id")
    assert len(count1) > 0
    assert count1[count1["quadkey"] == "12021023322022"]["count"].iloc[0] == 61
    assert count1[count1["quadkey"] == "12021023322023"]["count"].iloc[0] == 109
    assert count1[count1["quadkey"] == "12021023322032"]["count"].iloc[0] == 56


def test_get_adaptive_squares():
    # we use get_squares_polyfill. This method is used within aqk and provides information like children_id
    city_squares = get_squares_polyfill(city_polygon, 14)
    count1 = count_poi(city_squares, tess_data)

    assert len(get_adaptive_squares(count1, 100)) > len(city_squares)
    assert len(get_adaptive_squares(count1, len(tess_data))) == len(city_squares)


def test_split_linestring():
    road_data = pd.read_pickle("tests/road_data_Lille.pkl")
    split_road_data = split_linestring(road_data)

    assert len(split_road_data) > len(road_data)
