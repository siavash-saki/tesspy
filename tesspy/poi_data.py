import os
import numpy as np
import pandas as pd
import geopandas as gpd

import warnings
import osmnx as ox
import requests
import json
from shapely.geometry import Point


def geom_ceil(coordinate, precision=4):
    coordinate_ceil = np.true_divide(
        np.ceil(coordinate * 10**precision), 10**precision
    )
    return coordinate_ceil


def geom_floor(coordinate, precision=4):
    coordinate_floor = np.true_divide(
        np.floor(coordinate * 10**precision), 10**precision
    )
    return coordinate_floor


class POIdata:
    """
    This class creates a query for the investigated area and POI categories.
    The query is sent to osm using overpass API and the data is retrieved.

    Parameters
    ----------
    area : geopandas.GeoDataFrame or str
        GeoDataFrame must have a single shapely Polygon or MultiPolygon
        in geometry column and its CRS must be defined.
        str must be a name of a city, or an address of a region

    poi_categories : A list of OSM primary map features or 'all'
    timeout : int
        The TCP connection timeout for the overpass request
    verbose : bool
        If True, print information while computing
    """

    def __init__(self, area, poi_categories, timeout, verbose):
        self.area_buffered = None
        self.area = area
        self.poi_categories = poi_categories
        self.timeout = timeout
        self.verbose = verbose

    @staticmethod
    def osm_primary_features():
        """
        list of primary OSM features
        available at https://wiki.openstreetmap.org/wiki/Map_features

        Returns
        --------
        osm_primary_features_lst : list
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

    def create_overpass_query_string(self):
        """
        creates the query string to be passed to overpass

        Returns
        --------
        query_string : str
        """
        # make the query area a bit larger
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.area_buffered = self.area.buffer(0.008).simplify(0.005)

        exter_cooridinates = self.area_buffered.iloc[0].exterior.coords

        xy = np.array(exter_cooridinates)

        # make polygon string for OSM overpass query
        # Using the polygon, fewer data are retrieved, and it's faster but request is long can can lead to 414
        # poly_str = ''
        # for lat, lon in zip(xy[:, 1], xy[:, 0]):
        #     poly_str = poly_str + str(lat) + ' ' + str(lon) + ' '
        # poly_str = poly_str.strip()

        # make bounding box for OSM overpass query
        lat_min = geom_floor(np.min(xy[:, 0]))
        lon_min = geom_floor(np.min(xy[:, 1]))
        lat_max = geom_ceil(np.max(xy[:, 0]))
        lon_max = geom_ceil(np.max(xy[:, 1]))

        # if poi not in primary --> error
        # todo: is this necessary?
        for poi_category in self.poi_categories:
            if poi_category not in self.osm_primary_features():
                raise ValueError(
                    f"{poi_category} is not a valid POI primary category. "
                    f"See a list of OSM primary "
                    f"features with Tessellation.osm_primary_features()"
                )

        # create query string for overpass
        query_string = ""
        for element in ["node", "way"]:
            for poi_category in self.poi_categories:
                # query with polygon
                # query_string = query_string + f'{element}[{poi_category}](poly:"{poly_str}");'
                # query with bounding box
                query_string = query_string + f"{element}[{poi_category}];"

        # query_string = f"[out:json][timeout:{self.timeout}];(" + query_string + ');out geom;'
        query_string = (
            f"[bbox][out:json][timeout:{self.timeout}];("
            + query_string
            + ");out geom;"
            + f"&bbox={lat_min},{lon_min},{lat_max},{lon_max}"
        )

        return query_string

    def get_poi_data(self):
        """
        sends the query to osm using the overpass API and gets the data

        Returns
        --------
        poi_df : pandas.DataFrame
            A dataframe containing the POI, POI type, and coordinates
        """

        query_string = self.create_overpass_query_string()
        request_header = "https://overpass-api.de/api/interpreter?data="

        if self.verbose:
            print("Getting data from OSM...")

        # sending the request
        # todo: polygon instead bbox??
        # resp = requests.get("https://overpass-api.de/api/interpreter", data=query_string)
        resp = requests.get(url=request_header + query_string)
        if resp.status_code == 429:
            raise RuntimeError(
                "429 Too Many Requests:\n"
                "You have sent multiple requests from the same "
                "IP and passed the passed the fair use policy. "
                "Please wait a couple of minutes and then try again."
            )
        elif resp.status_code == 504:
            raise RuntimeError(
                "504 Gateway Timeout:\n"
                "the server has already so much load that"
                " the request cannot be executed."
                "Please try again later"
            )
        elif resp.status_code != 200:
            raise RuntimeError("Bad Request!")
        else:
            resp = json.loads(resp.text)

        if self.verbose:
            print("Creating POI DataFrame...")

        lst_nodes = []
        lst_ways = []

        generator = resp["elements"]
        for item in generator:

            for cat in self.poi_categories:
                if cat in item["tags"].keys():
                    item[cat] = True
            if item["type"] == "node":
                lst_nodes.append(item)
            elif item["type"] == "way":
                item["center_latitude"] = np.mean(
                    [point["lat"] for point in item["geometry"]]
                )
                item["center_longitude"] = np.mean(
                    [point["lon"] for point in item["geometry"]]
                )
                lst_ways.append(item)
            else:
                continue

        if self.verbose:
            print("Cleaning POI DataFrame...")

        nodes_df = pd.DataFrame(lst_nodes)
        ways_df = pd.DataFrame(lst_ways)

        if len(nodes_df) > 0 and len(ways_df) > 0:
            if self.verbose:
                print("Joining nodes and ways")

            nodes_df["geometry"] = nodes_df[["lon", "lat"]].apply(
                lambda p: [{"lat": p["lat"], "lon": p["lon"]}], axis=1
            )
            nodes_df = nodes_df.rename(
                columns={"lat": "center_latitude", "lon": "center_longitude"}
            )
            nodes_df = nodes_df.drop(columns=["id"])
            ways_df = ways_df.drop(columns=["id", "bounds", "nodes"])

            poi_df = pd.concat([ways_df, nodes_df]).fillna(False)

        elif len(nodes_df) == 0 and len(ways_df) > 0:
            if self.verbose:
                print("No nodes found. Return ways only.")

            ways_df = ways_df.drop(columns=["id", "bounds", "nodes"])
            poi_df = ways_df.fillna(False)

        elif len(nodes_df) > 0 and len(ways_df) == 0:
            if self.verbose:
                print("No ways found returning nodes only")

            nodes_df["geometry"] = nodes_df[["lon", "lat"]].apply(
                lambda p: [{"lat": p["lat"], "lon": p["lon"]}], axis=1
            )

            nodes_df = nodes_df.rename(
                columns={"lat": "center_latitude", "lon": "center_longitude"}
            )
            nodes_df = nodes_df.drop(columns=["id"])
            poi_df = nodes_df.fillna(False)
        else:
            raise ValueError(
                "No POI data for the inserted poi-category/categories and area."
            )

        # add POI categories with no data to poi_df
        for poi_category in self.poi_categories:
            if not hasattr(poi_df, poi_category):
                poi_df[poi_category] = False

        # prettify poi_df by sorting columns
        first_cols = ["type", "geometry", "tags", "center_latitude", "center_longitude"]

        second_cols = sorted(poi_df.columns.drop(first_cols))
        poi_df = poi_df[first_cols + second_cols]
        poi_df = poi_df.reset_index(drop=True)

        # Keep only the POI within the area
        geometry_column = [
            Point(coords)
            for coords in poi_df[["center_longitude", "center_latitude"]].values
        ]
        poi_geo_df = gpd.GeoDataFrame(geometry=geometry_column, crs="EPSG:4326")
        area_buffered_gdf = gpd.GeoDataFrame(
            geometry=self.area_buffered, crs="epsg:4326"
        )
        idx_to_keep = gpd.sjoin(poi_geo_df, area_buffered_gdf, predicate="within").index
        poi_df = poi_df.loc[idx_to_keep]

        if len(poi_df) == 0:
            raise ValueError("No poi data found that is within the area.")

        return poi_df


class RoadData:
    """
    This class creates a custom filter for the investigated area and highway types.
    The custom filter is used to collect road data using osmnx.

    Parameters
    ----------
    area : geopandas.GeoDataFrame or str
        GeoDataFrame must have a single shapely Polygon or MultiPolygon
        in geometry column and its CRS must be defined.
        str must be a name of a city, or an address of a region

    detail_deg : int
        integer to select the top (int) highway types
    split_roads : bool
        decide if LineStrings will be split up such that each LineString contains exactly 2 Points
    verbose : bool
        If True, print information while computing
    """

    def __init__(self, area, detail_deg=None, split_roads=True, verbose=False):

        self.detail_deg = detail_deg
        self.area = area
        self.verbose = verbose
        self.split_roads = split_roads

    @staticmethod
    def osm_highway_types():
        """
        list of OSM highway types
        available at https://wiki.openstreetmap.org/wiki/Key:highway

        Returns
        --------
        osm_highways_lst : list
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

    def create_custom_filter(self):
        """
        Creates custom filter that is used to collect road data with osmnx. The detail_deg is used to use only the top
        highway types.

        Returns
        --------
        custom_filter : list
            string with highway types that is used to collect data using osmnx
        """
        # Creating custom filter for all highway types
        if self.detail_deg is None:
            highwaytypes = self.osm_highway_types()

        # Creating custom filter for the specified detail degree
        elif type(self.detail_deg) is int:
            highwaytypes = self.osm_highway_types()[: self.detail_deg]

        else:
            raise ValueError("Please insert a valid detail degree: None or int")

        query = ""
        for types in highwaytypes[:-1]:
            query += f"{types}|"
        query += highwaytypes[-1]
        custom_filter = f"['highway'~'{query}']"

        if self.verbose:
            print(f"Selected highway type(s) are {custom_filter}")

        return custom_filter

    def get_road_network(self):
        """
        Collects the road network data based on the defined custom filter. The initial road data is based on graphs
        and transformed into a GeoDataFrame. Within the RoadData Class the user can divide to split the data such
        that each LineString (representing a street segment) contains only two points and not multiple points.

        Returns
        --------
        graph_edges_as_gdf : geopandas.GeoDataFrane
            GeoDataFrame containing road network
        """
        cf = self.create_custom_filter()
        if self.verbose:
            print("Collecting road network data...")
        graph = ox.graph_from_polygon(
            self.area.boundary.convex_hull.values[0], custom_filter=cf
        )
        graph_projected = ox.project_graph(graph, to_crs="epsg:4326")
        graph_undirected = graph_projected.to_undirected()
        graph_edges_as_gdf = ox.graph_to_gdfs(graph_undirected, nodes=False, edges=True)

        if self.verbose:
            print(f"Collected data has {len(graph_edges_as_gdf)} street segments.")

        return graph_edges_as_gdf
