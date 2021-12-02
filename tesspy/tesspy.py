import shapely
from shapely.geometry import Point, Polygon, LineString
import geopandas as gpd
import osmnx as ox
import h3pandas
from babelgrid import Babel


def add_numbers(a, b):
    return a + b


class TessObj():

    def __init__(self, city, resolution):
        self.resolution = resolution
        self.city = city
        pass

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
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :return:GeoDataFrame with all the squares
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = self.get_admin_Polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city
            pass
        else:
            print("Check available data formats")
        tiles = Babel('bing').polyfill(df_city.geometry.unary_union, resolution=resolution)


    def hexagon(self, city):
        pass

    def voronoi(self, city):
        pass

    def cityblocks(self, city):
        pass
