import os
import numpy as np
import pandas as pd
import geopandas as gpd

import warnings

from sklearn.cluster import KMeans
from scipy.spatial import Voronoi
import hdbscan

import tesspy
from tessellation_functions import *
from poi_data import *


# todo: return gdf with POI count


def get_city_polygon(city: str):
    """
    Gets the polygon of a city or an area

    Parameters
    ----------
    city : str
        city must be a name of a city, or an address of a region

    Returns
    --------
    df_city : geopandas.GeoDataFrame
        GeoDataFrame containing the polygon of the city
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df_city = ox.geocode_to_gdf(city)
    df_city = df_city[["osm_id", "geometry"]]
    df_city = df_city.rename({"osm_id": "osmid"})
    return df_city


def _check_input_geodataframe(gdf):
    """
    checks if the input gdf is in a correct format. Otherwise, an error is raised to
    demonstrate what the problem is.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        input GeoDataFrame

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        same GeoDataFrame as input
    """
    if len(gdf) != 1:
        raise ValueError("GeoDataFrame must have only one geometry element")
    else:
        if not hasattr(gdf, "geometry"):
            raise TypeError("Geometry column missing in GeoDataFrame")
        else:
            if type(gdf["geometry"].iloc[0]) not in [Polygon, MultiPolygon]:
                raise TypeError(
                    "Geometry must be of type shapely polygon or multipolygon"
                )
            else:
                if gdf.crs is None:
                    raise ValueError("GeoDataFrame must have a CRS")
                elif gdf.crs != "epsg:4326":
                    return gdf.to_crs(epsg=4326)
                else:
                    return gdf


def _check_valid_geometry_gdf(gdf):
    """
    Checks if the geometry type of the gdf is correct. We want only shapely.Polygon geom_types.
    If there are MultiPolygons in the input gdf, it will be exploded into Polygons

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame

    Returns
    --------
    gdf : geopandas.GeoDataFrame
    """

    if len(gdf) < 1:
        raise ValueError("GeoDataFrame must have only one geometry element")
    else:
        if not hasattr(gdf, "geometry"):
            raise TypeError("Geometry column missing in GeoDataFrame")
        else:
            print("MultiPolygon found. Splitting it up...")
            if "MultiPolygon" in gdf.geom_type.unique():
                gdf = explode(gdf)
                gdf = gdf.reset_index()
                gdf.drop(columns=["level_0", "level_1"], inplace=True)

                return gdf
            else:
                return gdf


def count_poi_per_tile(city, gdf, poi_categories=["amenity", "building"], timeout=120):
    """
    Counts different POI-categories per tile. For each POI-categories an additional column
    is added to the Tessellation GeoDataFrame (gdf). After counting each POI-category the
    final Tessellation GeoDataFrame with the additional count-columns is returned.

    Parameters
    ----------
    city: str
        Tessellation object or string of the underlying city
    gdf : geopandas.GeoDataFrame
        GeoDataFrame of the Tessellation for the underlying city
    poi_categories : list, default=["amenity", 'building']
        POI which will be count per tile. These POI may differ from the POI used
        to create the tessellation
        This can be a list of OSM primary map features or 'all'
            'all' means all the available POI categories
            Possible values in the list: ['aerialway', 'aeroway',
            'amenity', 'barrier', 'boundary', 'building', 'craft',
            'emergency', 'geological', 'healthcare', 'highway',
            'historic', 'landuse', 'leisure', 'man_made', 'military',
            'natural', 'office', 'place', 'power', 'public_transport',
            'railway', 'route', 'shop', 'sport', 'telecom', 'tourism',
            'water', 'waterway']
    timeout: int, default=120
        positive number indicating the time to wait for OSM to return POI data

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        Tessellation GeoDataFrame with additional columns (count_poi-columns are added)
    """

    if type(city) == str:
        city = Tessellation(city)

    else:
        raise ValueError(
            "Please insert a valid city-type. Valid types are: tesspy.Tessellation Obejct or String"
        )

    if len(gdf) < 1:
        raise ValueError(
            "Please insert a valid Tessellation-GeoDataFrame. A valid Tessellation-GeoDataFrame should "
            "include at least one tile."
        )
    if type(poi_categories) == str:
        poi_categories = [poi_categories]
    elif type(poi_categories) == list or type(poi_categories) == np.array():
        poi_categories = poi_categories
    else:
        raise ValueError(
            "Please insert valid poi_categories. Valid types are: string values based on "
            "osm_primary_features of list/numpy.array with osm_primary_features."
        )

    df_poi = POIdata(
        city.get_polygon(),
        poi_categories=poi_categories,
        timeout=timeout,
        verbose=False,
    ).get_poi_data()

    points_geom = df_poi[["center_longitude", "center_latitude"]].apply(
        lambda p: Point(p["center_longitude"], p["center_latitude"]), axis=1
    )

    tess_data = gpd.GeoDataFrame(
        geometry=points_geom, data=df_poi[poi_categories], crs="EPSG:4326"
    )

    tess_data["value"] = (
        tess_data.drop(columns=["geometry"]).idxmax(1).where(tess_data.any(1))
    )
    tess_data = tess_data[["value", "geometry"]]

    try:
        idx = [s for s in gdf.columns if s.__contains__("id")][0]
    except:
        idx = [s for s in gdf.columns if s.__contains__("key")][0]

    spatial_join = gpd.sjoin(gdf, tess_data)
    pivot_table = pd.pivot_table(
        spatial_join, index=idx, columns="value", aggfunc={"value": len}
    )

    pivot_table.columns = pivot_table.columns.droplevel()

    merged_polygons = gdf.merge(pivot_table, how="left", on=idx)
    merged_polygons.fillna(0, inplace=True)

    return merged_polygons


class Tessellation:
    """
    Creates a Tessellation object using a GeoDataFrame or a city name,
    which allows for different tessellation methods

    Parameters
    ----------
    area : geopandas.GeoDataFrame or str
        GeoDataFrame must have a single shaply Polygon or MultiPolygon
        in geometry column and its CRS must be defined.
        str must be a name of a city, or an address of a region

    Examples
    --------
    >>> ffm= Tessellation('Frankfurt am Main')
    """

    def __init__(self, area):

        if type(area) == gpd.GeoDataFrame:
            self.area_gdf = _check_input_geodataframe(area)
        elif type(area) == str:
            self.area_gdf = get_city_polygon(area)
        else:
            raise TypeError("area must be in format: GeoDataFrame or string")

        self.poi_dataframe = pd.DataFrame()
        self.available_poi_categories = []
        self.road_network = pd.DataFrame()
        self.queried_highway_types = []

    def _get_missing_poi_categories(self, poi_list):
        """
        Checks if the poi categories are already available in the object poi_dataframe
        Creates a list of missing categories, which should be downloaded from OSM

        Parameters
        ----------
        poi_list : list
            A list of passed POI categories

        Returns
        -------
        missing_poi_categories : list
            A list containing the POI categories which should be downloaded
        """

        missing_poi_categories = []
        for poi_category in poi_list:
            if not hasattr(self.poi_dataframe, poi_category):
                missing_poi_categories.append(poi_category)
        return missing_poi_categories

    def squares(self, resolution: int):
        """
        Generate square grid laying over the area

        Parameters
        ----------
        resolution : int
            Specifies the size of squares
            A positive number between

        Returns
        -------
        df_qk_squares : pandas.DataFrame
            Dataframe containing squares
        """
        df_qk_squares = get_squares_polyfill(self.area_gdf, resolution)
        df_qk_squares = df_qk_squares.drop(columns=["osm_id", "children_id"])

        return df_qk_squares

    def hexagons(self, resolution: int):
        """
        Generate hexagon grid laying over the area

        Parameters
        ----------
        resolution : int
            Specifies the size of hexagons
            A positive number between

        Returns
        -------
        df_h3_hexagons : pandas.DataFrame
            Dataframe containing hexagons
        """

        df_h3_hexagons = get_h3_hexagons(self.area_gdf, resolution)
        df_h3_hexagons = df_h3_hexagons.reset_index().rename(
            columns={"index": "hex_id"}
        )

        return df_h3_hexagons

    def adaptive_squares(
        self,
        start_resolution: int,
        poi_categories=["amenity", "building"],
        threshold=None,
        timeout=60,
        verbose=False,
    ):

        """
        Generate adaptive squares based on the input POI data.
        Squares are created at the start resolution. Each square is broken
        into four smaller squares while the number of its POI exceeds the
        threshold.
        POI categories should be a list of OSM primary map features.
        A complete list can be found on the OSM website:
        https://wiki.openstreetmap.org/wiki/Map_features

        Parameters
        ----------
        start_resolution : int
            Specifies the size of initial squares
        poi_categories : A list of OSM primary map features or 'all'
                         default=["amenity", 'building']
            'all' means all the available POI categories
            Possible values in the list: ['aerialway', 'aeroway',
            'amenity', 'barrier', 'boundary', 'building', 'craft',
            'emergency', 'geological', 'healthcare', 'highway',
            'historic', 'landuse', 'leisure', 'man_made', 'military',
            'natural', 'office', 'place', 'power', 'public_transport',
            'railway', 'route', 'shop', 'sport', 'telecom', 'tourism',
            'water', 'waterway']
        threshold : int, default=None
            Threshold for the number of POI in a single square. If square
            has more, it is divided into four squares. If None passed, the
            median number of POI per square in the initial level is used
            as threshold.
        timeout : int, default=60
            The TCP connection timeout for the request
        verbose : bool, default=False
            If True, print information while computing

        Returns
        -------
        df_adaptive_squares : pandas.DataFrame
            Dataframe containing adaptive squares
        """

        if poi_categories == "all":
            poi_categories = self.osm_primary_features()

        # check if the categories are already available in poi_data
        # find the ones which are not available and should be downloaded from OSM
        missing_poi_categories = self._get_missing_poi_categories(poi_categories)

        # if there is any category that should be downloaded -->
        # create the query and run it
        if len(missing_poi_categories) > 0:
            poi_data_obj = POIdata(
                self.area_gdf, missing_poi_categories, timeout, verbose
            )
            poi_data_new = poi_data_obj.get_poi_data()
            # concat the data with the available poi_data dataframe of the Tessellation object
            self.poi_dataframe = pd.concat([self.poi_dataframe, poi_data_new]).fillna(
                False
            )
            self.poi_dataframe = self.poi_dataframe.reset_index(drop=True)

        # data, based on which, tessellation should be done
        tess_data = self.poi_dataframe[
            self.poi_dataframe[poi_categories].sum(axis=1) > 0
        ]
        points_geom = tess_data[["center_longitude", "center_latitude"]].apply(
            lambda p: Point(p["center_longitude"], p["center_latitude"]), axis=1
        )
        tess_data = gpd.GeoDataFrame(
            geometry=points_geom, data=tess_data[poi_categories], crs="EPSG:4326"
        )
        poi_data_aqk = tess_data.rename(columns={"points_geom": "geometry"})
        df_aqk = get_squares_polyfill(self.area_gdf, start_resolution)
        aqk_count_df = count_poi(df_aqk, poi_data_aqk)

        if not threshold:
            threshold = int(np.median(aqk_count_df["count"].values))
            if verbose:
                print(
                    f"Threshold={threshold}  ==> set as the median POI-count per square at the initial level"
                )

        i = start_resolution
        while max(aqk_count_df["count"].values) > threshold:
            i += 1
            if verbose:
                print(f"Threshold exceeded! Squares are subdivided into resolution {i}")

            df_tmp = get_adaptive_squares(aqk_count_df, threshold)
            df_tmp.drop(columns=["count"], inplace=True)
            df_tmp2 = count_poi(df_tmp, poi_data_aqk)
            aqk_count_df = df_tmp2

        final_aqk = gpd.sjoin(aqk_count_df, self.area_gdf)
        final_aqk = final_aqk.drop(columns=["osm_id", "children_id", "index_right"])

        return final_aqk

    def voronoi(
        self,
        cluster_algo="k-means",
        poi_categories=["amenity", "building"],
        timeout=60,
        n_polygons=100,
        min_cluster_size=15,
        verbose=False,
    ):
        """
        Generate Voronoi polygons based on the input POI data.
        POI categories should be a list of OSM primary map features.
        A complete list can be found on the OSM website:
        https://wiki.openstreetmap.org/wiki/Map_features

        Parameters
        ----------
        cluster_algo : {'k-means', 'hdbscan', None}, default='k-means'
            Algorithm for clustering the POI data before creating
            Voronoi generators. If None passed, POI data are
            directly used as generators.
        poi_categories : A list of OSM primary map features or 'all'
                         default=["amenity", 'building']
            'all' means all the available POI categories
            Possible values in the list: ['aerialway', 'aeroway',
            'amenity', 'barrier', 'boundary', 'building', 'craft',
            'emergency', 'geological', 'healthcare', 'highway',
            'historic', 'landuse', 'leisure', 'man_made', 'military',
            'natural', 'office', 'place', 'power', 'public_transport',
            'railway', 'route', 'shop', 'sport', 'telecom', 'tourism',
            'water', 'waterway']
        timeout : int, default=60
            The TCP connection timeout for the request
        n_polygons : int, default=100
            Only when cluster_algo="k-means", approximate number of
            polygons to be created. Positive number and less than
            the initial POI numbers
        min_cluster_size : int, default=15
            Only when cluster_algo="hdbscan", minimum cluster size
            Positive number
        verbose : bool, default=False
            If True, print information while computing

        Returns
        -------
        df_voronoi : pandas.DataFrame
            Dataframe containing Voronoi polygons
        """

        if poi_categories == "all":
            poi_categories = self.osm_primary_features()

        # check if the categories are already available in poi_data
        # find the ones which are not available and should be downloaded from OSM
        missing_poi_categories = self._get_missing_poi_categories(poi_categories)

        if type(self.area_gdf) == MultiPolygon:
            queried_area = self.area_gdf.convex_hull
        else:
            queried_area = self.area_gdf

        # if there is any category that should be downloaded -->
        # create the query and run it
        if len(missing_poi_categories) > 0:
            poi_data_obj = POIdata(
                queried_area, missing_poi_categories, timeout, verbose
            )
            poi_data_new = poi_data_obj.get_poi_data()
            # concat the data with the available poi_data dataframe of the Tessellation object
            self.poi_dataframe = pd.concat([self.poi_dataframe, poi_data_new]).fillna(
                False
            )
            self.poi_dataframe = self.poi_dataframe.reset_index(drop=True)

        # data, based on which, tessellation should be done
        tess_data = self.poi_dataframe[
            self.poi_dataframe[poi_categories].sum(axis=1) > 0
        ]
        data_locs = tess_data[["center_longitude", "center_latitude"]].values

        # create generators for Voronoi diagram
        if cluster_algo == "k-means":
            if verbose:
                print("K-Means Clustering...")
            clustering = KMeans(n_clusters=n_polygons).fit(data_locs)
            generators = [
                np.mean(data_locs[clustering.labels_ == label], axis=0)
                for label in range(n_polygons)
            ]

        elif cluster_algo == "hdbscan":
            if verbose:
                print("HDBSCAN Clustering... This can take a while...")
            clustering = hdbscan.HDBSCAN(
                min_cluster_size=min_cluster_size, prediction_data=True
            ).fit(data_locs)
            generators = [
                np.mean(data_locs[clustering.labels_ == label], axis=0)
                for label in range(clustering.labels_.max() + 1)
            ]

        elif cluster_algo is None:
            if len(tess_data) > 5000:
                raise ValueError(
                    "Too many generators for Voronoi diagram. Please select a clustering algorithm"
                )
            else:
                generators = data_locs

        else:
            raise ValueError(
                "Please use a clustering algorithm: k-means, hdbscan or None"
            )

        # create Voronoi polygons
        if verbose:
            print("Creating Voronoi polygons...")
        voronoi_dia = Voronoi(generators)
        voronoi_poly = gpd.GeoDataFrame(
            geometry=[p for p in voronoi_polygons(voronoi_dia, 0.1)], crs="EPSG:4326"
        )
        voronoi_poly = gpd.sjoin(voronoi_poly, self.area_gdf)
        vor_polygons = voronoi_poly.intersection(self.area_gdf.geometry.iloc[0])
        df_voronoi = gpd.GeoDataFrame(geometry=vor_polygons)

        df_voronoi = _check_valid_geometry_gdf(df_voronoi)

        df_voronoi.reset_index(inplace=True)
        df_voronoi.rename(columns={"index": "voronoi_id"}, inplace=True)
        df_voronoi["voronoi_id"] = "voronoiID" + df_voronoi["voronoi_id"].astype(str)

        return df_voronoi

    def city_blocks(
        self, n_polygons=None, detail_deg=None, split_roads=True, verbose=False
    ):
        """
        Create city bocks (tiles) using road data from the area.
        To collect road data, specify the highway types by
        modifying detail_deg

        Parameters
        ----------
        n_polygons: int, default = None
            targeted number of city blocks, this is an approximation
            the final number of polygons can vary slightly
        detail_deg: int, default = None
            define the number of the top (osm) highway types to use for creating city blocks
        split_roads: bool, default = True
            if True, LineStrings are split up such that each LineString contains exactly 2 Points
            This usually make the polygonizing more robust but slower.
        verbose : bool, default = False
            If True, print information while computing

        Returns
        -------
        final_city_blocks : geopandas.GeoDataFrame
            GeoDataFrame with city block tiles
        """

        # check if the road data is already available
        if detail_deg is None:
            highwaytypes = self.osm_highway_types()
        elif type(detail_deg) is int and detail_deg <= len(self.osm_highway_types()):
            highwaytypes = self.osm_highway_types()[:detail_deg]
        else:
            raise ValueError("Please insert a valid detail degree: None or int")

        # for OSM query, we need a polygon. If multipolygon is passed, we generate the convex hull
        if type(self.area_gdf) == MultiPolygon:
            queried_area = self.area_gdf.convex_hull
        else:
            queried_area = self.area_gdf

        # if the road network is not available
        if self.queried_highway_types != highwaytypes:
            # collect road network data
            road_data_collect_object = RoadData(
                queried_area, detail_deg, split_roads, verbose
            )
            road_data = road_data_collect_object.get_road_network()
            self.road_network = road_data
            # keep track of downloaded road network to prevent similar OSM request
            self.queried_highway_types = highwaytypes
        # if road network available, don't download it again.
        else:
            road_data = self.road_network

        # todo: causing problems (kernel dies): should this stay?
        # split roads if True
        if split_roads:
            if verbose:
                print(
                    "Splitting the linestring, such that each linestring has exactly 2 points."
                )
            road_data = split_linestring(road_data)

        if verbose:
            print("Creating initial city blocks using the road network data...")

        # creating blocks
        blocks = create_blocks(road_data)

        # keep polygons inside studied area
        polygons_in_area = gpd.sjoin(blocks, queried_area, how="inner")
        polygons_in_area.drop(columns=["index_right"], inplace=True)

        # create polygons by the border
        rest_polygons = get_rest_polygon(polygons_in_area, queried_area)

        # add rest polygons to all polygons
        city_blocks = pd.concat([polygons_in_area, rest_polygons])

        # merging small polygons using hierarchical clustering
        if not n_polygons:
            city_blocks = city_blocks[["geometry"]].reset_index(drop=True)

            city_blocks = _check_valid_geometry_gdf(city_blocks)

            city_blocks.reset_index(inplace=True)
            city_blocks.rename(columns={"index": "cityblock_id"}, inplace=True)
            city_blocks["cityblock_id"] = "cityblockID" + city_blocks[
                "cityblock_id"
            ].astype(str)

            return city_blocks

        if n_polygons > len(city_blocks):
            raise ValueError(
                f"Cannot extract more city blocks than initial city blocks!"
                f"Initial city blocks are {len(city_blocks)}, desired city blocks are {n_polygons}"
                f"Choose a value less than {len(city_blocks)}."
            )

        if verbose:
            print("Merging small city blocks...")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            city_blocks["centroid"] = city_blocks.centroid

        coordinates = np.column_stack(
            [city_blocks["centroid"].x, city_blocks["centroid"].y]
        )
        # the algorithm needs O(n²) memory and O(n³) runtime
        # => doesn't work with large data ==> not enough RAM ==> kernel dies
        # think about changing the clustering algorithm
        model = AgglomerativeClustering(n_clusters=n_polygons, affinity="euclidean")
        model.fit(coordinates)

        city_blocks["Cluster"] = model.labels_

        merged_polys = []
        for idx in city_blocks["Cluster"].unique():
            tmp = city_blocks[city_blocks["Cluster"] == idx]
            polygons = tmp["geometry"].to_numpy()
            merged_polygon = gpd.GeoSeries(unary_union(polygons))
            merged_polys.append(merged_polygon[0])

        new_merged_polys = {"geometry": merged_polys}
        merged_polys_df = gpd.GeoDataFrame(new_merged_polys, crs="EPSG:4326")
        keep_df = merged_polys_df[merged_polys_df.geom_type == "Polygon"]
        to_explode = merged_polys_df[merged_polys_df.geom_type == "MultiPolygon"]
        explode_df = explode(to_explode)
        explode_df = explode_df.reset_index()
        explode_df.drop(columns=["level_0", "level_1"], inplace=True)
        final_city_blocks = pd.concat([keep_df, explode_df])
        final_city_blocks = final_city_blocks[["geometry"]].reset_index(drop=True)

        final_city_blocks = _check_valid_geometry_gdf(final_city_blocks)

        final_city_blocks.reset_index(inplace=True)
        final_city_blocks.rename(columns={"index": "cityblock_id"}, inplace=True)
        final_city_blocks["cityblock_id"] = "cityblockID" + final_city_blocks[
            "cityblock_id"
        ].astype(str)

        return final_city_blocks

    def get_polygon(self):
        """
        Returns
        --------
        area_gdf : geopandas.GeoDataFrame
            the area polygon in GeoDataFrame format
        """
        return self.area_gdf

    def get_poi_data(self):
        """
        Returns
        --------
        area_gdf : geopandas.GeoDataFrame
            the POI data in GeoDataFrame format
        """
        return self.poi_dataframe

    def get_road_network(self):
        """
        Returns
        --------
        road_network : geopandas.GeoDataFrame
            the road network data in GeoDataFrame format
        """
        return self.road_network

    @staticmethod
    def osm_primary_features():
        """
        list of primary OSM features
        available at https://wiki.openstreetmap.org/wiki/Map_features

        Returns
        --------
        osm_primary_features_lst: list
        """
        osm_primary_features_lst = [
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

        return osm_primary_features_lst

    @staticmethod
    def osm_highway_types():
        """
        list of all highway types

        Returns
        --------
        osm_highways_lst: list
            list of highway types
        """
        osm_highways_lst = [
            "motorway",
            "trunk",
            "primary",
            "secondary",
            "tertiary",
            "residential",
            "unclassified",
            "motorway_link",
            "trunk_link",
            "primary_link",
            "secondary_link",
            "living_street",
            "pedestrian",
            "track",
            "bus_guideway",
            "footway",
            "path",
            "service",
            "cycleway",
        ]

        return osm_highways_lst
