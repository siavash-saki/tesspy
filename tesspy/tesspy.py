import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox
# import h3
from babelgrid import Babel
import shapely
from shapely.geometry import Point, Polygon, LineString


def get_admin_polygon(city: str):
    """
    :param city:
    :return: simple GeoDataFrame containing only the Polygon of the city
    """
    try:
        df_city = ox.geocode_to_gdf(city)
        df_city = df_city[["osm_id", "geometry"]]
        df_city = df_city.rename({"osm_id": "osmid"})
        return df_city

    except:
        print("Input must be a city in string format")


def get_bbox(city: str):
    """

    :param city:
    :return: Array of Int which represent the bounding boy for the given city in the format
    [north, south,east,west]
    """
    gdf_city = ox.geocode_to_gdf(city)
    bbox = [gdf_city["bbox_north"].values[0],
            gdf_city["bbox_south"].values[0],
            gdf_city["bbox_east"].values[0],
            gdf_city["bbox_west"].values[0]]
    return bbox


def split_linestring(df):
    """

    :param df:
    :return: DataFrame with shapely.linestring with two points each
    The linestring, which are split up will have the same osmid since its still the same street
    """
    linestrings = []
    osmid = []

    for idx, row in df.iterrows():
        if len(row["geometry"].coords) == 2:
            linestrings.append(row["geometry"])
            osmid.append(row["osmid"])
        else:
            for i in range(0, len(row["geometry"].coords) - 1):
                p1 = Point(row["geometry"].coords[i][0], row["geometry"].coords[i][1])
                p2 = Point(row["geometry"].coords[i + 1][0], row["geometry"].coords[i + 1][1])
                linestrings.append(LineString([p1, p2]))
                osmid.append(row["osmid"])

    dataset = gpd.GeoDataFrame({"osmid": osmid, "geometry": linestrings})
    return dataset


def explode(gdf):
    """

    :param gdf: which can have multi-part geometries that will be exploded
    :return: GeoDataFrame with single geometries
    example: Multipolygon -> multiple Polygons
    """
    gs = gdf.explode()
    gdf2 = gs.reset_index().rename(columns={0: "geometry"})
    gdf_out = gdf2.merge(gdf.drop("geometry", axis=1), left_on="level_0", right_index=True)
    gdf_out = gdf_out.set_index(["level_0", "level_1"]).set_geometry("geometry")
    gdf_out.crs = gdf.crs
    return gdf_out


def count_poi(df, points):
    """
    :param df: Geodataframe of LGUs
    :param points: GeoDataFrame of POI
    :return: GeoDataFrame of LGUs with "count" column that counts points per LGU
    """

    pointsInPolygon = gpd.sjoin(df, points, how="left", op='contains')
    pointsInPolygon['count'] = 1
    pointsInPolygon.reset_index(inplace=True)
    tmp_a = pointsInPolygon.groupby(by='tile_id').count()['count'].reset_index().sort_values(by="tile_id",
                                                                                             ascending=True)
    tmp_b = df.reset_index().sort_values(by="tile_id", ascending=True)
    final_df = pd.merge(tmp_a, tmp_b, on="tile_id")
    final_df.set_index("tile_id", inplace=True)

    return final_df


def adaptive_tessellation(gdf, threshold):
    """
    :param gdf: GeoDataFrame from babel (with children_id column)
    :param threshold: Threshold of max data points per LGU allowed
    :return: divided GeoDataFrame where LGUs are subdivided that exceed the threshold
    """

    # DataFrame where the squares will be subdivided
    gdf_exceeded = gdf[gdf["count"] >= threshold]

    for idx in gdf_exceeded.index:
        # get all children (so array of length 4)
        children = gdf_exceeded.loc[[str(idx)]]["children_id"].values[0]

        # Delete the row in the overall dataframe
        gdf.drop([idx], inplace=True)

        for child in children:
            child_dict = Babel('bing').id_to_tile(child).to_dict()
            child_gdf = gpd.GeoDataFrame.from_dict(child_dict, orient='index')
            child_gdf = child_gdf.transpose()
            child_gdf = child_gdf.set_geometry("shapely")
            child_gdf = child_gdf.set_crs("EPSG:4326")
            child_gdf.set_index('tile_id', inplace=True)
            gdf = gdf.append(child_gdf)

    return gdf


class TessObj:

    def __init__(self, city, resolution):
        self.resolution = resolution
        self.city = city
        pass

    def quadKey(self, city, resolution: int):
        """
        :param resolution: resolution for mircosoft bing. Int ∈ [1,...,23]
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :return:GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")
        tiles = Babel('bing').polyfill(df_city.geometry.unary_union, resolution=resolution)
        df_qk = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely', crs="EPSG:4326")
        return df_qk

    def hexagon(self, city, resolution: int):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param resolution: resolution for mircosoft bing. Int ∈ [0,...,15]
        :return: GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")
        tiles = Babel('h3').polyfill(df_city.geometry.unary_union, resolution=resolution)
        df_h3 = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely', crs="EPSG:4326")
        return df_h3

    def adaptive_quadkey(self, city, poi_data, start_resolution):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param poi_data: GeoDataFrame of Data that will be used as subdivide criteria
        :param start_resolution: start resolution from which the algorithm will subdivide the area
        :return: GeoDataFrame of the LGUs with additional column "count" that counts the poi_data per LGU
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")

        tiles = Babel("bing").polyfill(df_city.geometry.unary_union, resolution=start_resolution)
        df_aqk = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely')
        df_aqk = df_aqk.set_crs("EPSG:4326")
        df_aqk.set_index("tile_id", inplace=True)

        aqk_count = count_poi(df_aqk, poi_data)
        aqk_count = gpd.GeoDataFrame(aqk_count, geometry="shapely")
        threshold = int(np.median(aqk_count["count"].values))
        while max(aqk_count["count"].values) > threshold:
            df_temp = adaptive_tessellation(aqk_count, threshold)
            df_temp.drop(columns=["count"], inplace=True)

            df_temp2 = count_poi(df_temp, poi_data)
            aqk_count = gpd.GeoDataFrame(df_temp2, geometry="shapely")

        return aqk_count

    def voronoi(self, city):
        pass

    def cityblocks(self, city):
        pass
