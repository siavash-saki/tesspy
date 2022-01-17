import hdbscan
import numpy as np
import pandas as pd
import geopandas as gpd
from geopandas import sjoin
import osmnx as ox
import h3
from babelgrid import Babel
import shapely
from shapely.geometry import Point, Polygon, LineString, mapping, MultiPoint
from shapely.ops import polygonize, cascaded_union
from sklearn.cluster import AgglomerativeClustering, KMeans
from scipy.spatial import Voronoi
from collections import defaultdict
import mercantile


def get_city_polygon(city: str):
    """
    :param city:
    :return: simple GeoDataFrame containing only the Polygon of the city
    """
    df_city = ox.geocode_to_gdf(city)
    df_city = df_city[["osm_id", "geometry"]]
    df_city = df_city.rename({"osm_id": "osmid"})
    return df_city



## we don't need this --> input either geopandas or string
# def get_bbox(city: str):
#     """
#     :param city:
#     :return: Array of Int which represent the bounding boy for the given city in the format
#     [north, south,east,west]
#     """
#     gdf_city = ox.geocode_to_gdf(city)
#     bbox = [gdf_city["bbox_north"].values[0],
#             gdf_city["bbox_south"].values[0],
#             gdf_city["bbox_east"].values[0],
#             gdf_city["bbox_west"].values[0]]
#     return bbox
#
#
#     dataset = gpd.GeoDataFrame({"osmid": osmid, "geometry": linestrings})
#     return dataset




class Tessellation:

    def __init__(self, area):

        # TODO: are should be a geodataframe with a single polygone (or MultiPolygon) and available crs
        # TODO: if crs is not epsg:4326 then convert
        if type(area) == gpd.GeoDataFrame:
            self.area_gdf = area

        elif type(area) == str:
            self.area_gdf = get_city_polygon(area)

        else:
            raise TypeError("area must be in format: GeoDataFrame or string")

        self.poi_data = None
        self.road_network = None

    def squares(self, resolution: int):
        """
        :param resolution: resolution for mircosoft bing. Int ∈ [1,...,23]
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :return:GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """
        df_qk = square_polyfill(self.area_gdf, resolution)
        return df_qk

    def hexagons(self, resolution: int):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param resolution: resolution for mircosoft bing. Int ∈ [0,...,15]
        :return: GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """

        df_h3 = get_h3(self.area_gdf, resolution)
        return df_h3

    def adaptive_quadkey(self, start_resolution, poi_categories=["amenity",'building'], threshold=None):


        # tiles = Babel("bing").polyfill(df_city.geometry.unary_union, resolution=start_resolution)
        # df_aqk = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely')
        # df_aqk = df_aqk.set_crs("EPSG:4326")
        # df_aqk.set_index("tile_id", inplace=True)
        #
        # aqk_count = count_poi(df_aqk, poi_data)
        # aqk_count = gpd.GeoDataFrame(aqk_count, geometry="shapely")
        # threshold = int(np.median(aqk_count["count"].values))
        # while max(aqk_count["count"].values) > threshold:
        #     df_temp = adaptive_tessellation(aqk_count, threshold)
        #     df_temp.drop(columns=["count"], inplace=True)
        #
        #     df_temp2 = count_poi(df_temp, poi_data)
        #     aqk_count = gpd.GeoDataFrame(df_temp2, geometry="shapely")

        return aqk_count

    def voronoi(self, cluster_algo="k-means"):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param data: GeoDataFrame with a geometry column representing data like POI or generators for Voronoi Diagrams
        :param cluster_algo: if not generators are used directly the spatial data is going to be clustered. specify the
        cluster algorithm
        :return: GeoDataFrame with Polygons representing the voronoi polygons
        """

        data_locs = np.column_stack([data.geometry.x, data.geometry.y])

        # TODO: maybe include more cluster algorithms like DBSCAN,..
        if cluster_algo == "k-means":
            # TODO: Think about automating nb of clusters
            clusterer = KMeans(n_clusters=500).fit(data_locs)
            data["ClusterIDs"] = clusterer.labels_

        elif cluster_algo == "hdbscan":
            # TODO: Think about automating min_cluster_size e.g. x% of the initial data shape
            clusterer = hdbscan.HDBSCAN(min_cluster_size=15, prediction_data=True).fit(data_locs)
            soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
            data["ClusterIDs"] = [np.argmax(x) for x in soft_clusters]

        elif cluster_algo == "None":
            data["ClusterIDs"] = data.geometry

        else:
            raise ValueError("Please use one of the implemented clustering: k-mean, hdbscan, dbscan")

        convex_hulls = [MultiPoint(data[data["ClusterIDs"] == ID]["geometry"].values).convex_hull for ID in
                        data["ClusterIDs"].unique()]
        gdf_generators = gpd.GeoDataFrame(geometry=convex_hulls, crs="EPSG:4326")
        generators = gdf_generators["geometry"].centroid

        voronoi_dia = Voronoi(np.column_stack([generators.x, generators.y]))
        voronoi_poly = gpd.GeoDataFrame(geometry=[p for p in voronoi_polygons(voronoi_dia, 0.1)], crs="EPSG:4326")
        vor_polygons = voronoi_poly.intersection(df_city.geometry.iloc[0])

        df_vor = gpd.GeoDataFrame(geometry=vor_polygons)
        return df_vor

    def city_blocks(self):
        pass
