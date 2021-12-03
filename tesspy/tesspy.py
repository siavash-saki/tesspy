
import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
# import h3
from babelgrid import Babel
import shapely
from shapely.geometry import Point, Polygon, LineString


def add_numbers(a, b):
    return a + b


class TessObj:

    def __init__(self, city, resolution):
        self.resolution = resolution
        self.city = city
        pass

#should be define this function external s.t. its not a class object?
    def get_admin_Polygon(self, city):
        try:
            df_city = ox.geocode_to_gdf(city)
            df_city = df_city[["osm_id", "geometry"]]
            df_city = df_city.rename({"osm_id": "osmid"})
            return df_city

        except:
            print("Input must be a city in string format")

    def quadKey(self, city, resolution):
        """
        :param resolution: resolution for mircosoft bing. Int ∈ [1,...,23]
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :return:GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = self.get_admin_Polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")
        tiles = Babel('bing').polyfill(df_city.geometry.unary_union, resolution=resolution)
        df_qk = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely', crs="EPSG:4326")
        return df_qk

    def hexagon(self, city, resolution):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param resolution: resolution for mircosoft bing. Int ∈ [0,...,15]
        :return: GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = self.get_admin_Polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")
        tiles = Babel('bing').polyfill(df_city.geometry.unary_union, resolution=resolution)
        df_h3 = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely', crs="EPSG:4326")
        return df_h3

    def voronoi(self, city):
        pass

    def cityblocks(self, city):
        pass
