"""
POI (Point of Interest) data retrieval via the OSM Overpass API.
"""

import json
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from shapely.geometry import Point

from tesspy._constants import OSM_PRIMARY_FEATURES
from tesspy.data._overpass import geom_ceil, geom_floor


class POIdata:
    """
    Query the OSM Overpass API for Points of Interest within a study area.

    Parameters
    ----------
    area : geopandas.GeoDataFrame
        GeoDataFrame with a single Polygon or MultiPolygon and a defined CRS.
    poi_categories : list of str
        OSM primary map feature categories to query.
    timeout : int
        TCP connection timeout in seconds for the Overpass request.
    verbose : bool
        If True, print progress information.
    """

    def __init__(
        self,
        area: gpd.GeoDataFrame,
        poi_categories: list[str],
        timeout: int,
        verbose: bool,
    ) -> None:
        self.area_buffered = None
        self.area = area
        self.poi_categories = poi_categories
        self.timeout = timeout
        self.verbose = verbose

    @staticmethod
    def osm_primary_features() -> list[str]:
        """
        Return the list of primary OSM map feature categories.
        See https://wiki.openstreetmap.org/wiki/Map_features

        Returns
        --------
        list of str
        """
        return OSM_PRIMARY_FEATURES

    def create_overpass_query_string(self) -> str:
        """
        Build the Overpass API query string for the study area.

        Returns
        --------
        query_string : str
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.area_buffered = self.area.buffer(0.008).simplify(0.005)

        exter_coordinates = self.area_buffered.iloc[0].exterior.coords
        xy = np.array(exter_coordinates)

        lat_min = geom_floor(np.min(xy[:, 0]))
        lon_min = geom_floor(np.min(xy[:, 1]))
        lat_max = geom_ceil(np.max(xy[:, 0]))
        lon_max = geom_ceil(np.max(xy[:, 1]))

        for poi_category in self.poi_categories:
            if poi_category not in self.osm_primary_features():
                raise ValueError(
                    f"{poi_category} is not a valid POI primary category. "
                    f"See a list of OSM primary features with "
                    f"Tessellation.osm_primary_features()"
                )

        query_string = ""
        for element in ["node", "way"]:
            for poi_category in self.poi_categories:
                query_string = query_string + f"{element}[{poi_category}];"

        query_string = (
            f"[bbox][out:json][timeout:{self.timeout}];("
            + query_string
            + ");out geom;"
            + f"&bbox={lat_min},{lon_min},{lat_max},{lon_max}"
        )

        return query_string

    def get_poi_data(self) -> pd.DataFrame:
        """
        Send the Overpass query and parse the returned POI data.

        Returns
        --------
        poi_df : pandas.DataFrame
            DataFrame with POI type, geometry, tags, center coordinates,
            and boolean columns for each queried POI category.
        """
        query_string = self.create_overpass_query_string()
        request_header = "https://overpass-api.de/api/interpreter?data="

        if self.verbose:
            print("Getting data from OSM...")

        resp = requests.get(url=request_header + query_string)
        if resp.status_code == 429:
            raise RuntimeError(
                "429 Too Many Requests:\n"
                "You have sent multiple requests from the same IP and exceeded "
                "the fair use policy. Please wait a few minutes and try again."
            )
        elif resp.status_code == 504:
            raise RuntimeError(
                "504 Gateway Timeout:\n"
                "The server is under heavy load and cannot process the request. "
                "Please try again later."
            )
        elif resp.status_code != 200:
            raise RuntimeError("Bad Request!")
        else:
            resp = json.loads(resp.text)

        if self.verbose:
            print("Creating POI DataFrame...")

        lst_nodes = []
        lst_ways = []

        for item in resp["elements"]:
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
                print("No nodes found. Returning ways only.")

            ways_df = ways_df.drop(columns=["id", "bounds", "nodes"])
            poi_df = ways_df.fillna(False)

        elif len(nodes_df) > 0 and len(ways_df) == 0:
            if self.verbose:
                print("No ways found. Returning nodes only.")

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
                "No POI data found for the specified poi_categories and area."
            )

        for poi_category in self.poi_categories:
            if not hasattr(poi_df, poi_category):
                poi_df[poi_category] = False

        first_cols = ["type", "geometry", "tags", "center_latitude", "center_longitude"]
        second_cols = sorted(poi_df.columns.drop(first_cols))
        poi_df = poi_df[first_cols + second_cols]
        poi_df = poi_df.reset_index(drop=True)

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
            raise ValueError("No POI data found within the study area.")

        return poi_df
