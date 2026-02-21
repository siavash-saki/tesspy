"""
POI (Point of Interest) data retrieval via the OSM Overpass API.
"""

import json
import logging
import time
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from shapely.geometry import Point

from tesspy._constants import OSM_PRIMARY_FEATURES
from tesspy._logging import log_progress
from tesspy.data._overpass import geom_ceil, geom_floor

logger = logging.getLogger(__name__)


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
        If True, log progress information via the ``tesspy`` logger.
    """

    def __init__(
        self,
        area: gpd.GeoDataFrame,
        poi_categories: list[str],
        timeout: int,
        verbose: bool,
    ) -> None:
        self.area_buffered: gpd.GeoSeries | None = None
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
        query_build_start = time.perf_counter()
        log_progress(
            logger,
            self.verbose,
            "event=poi.query.build.start poi_categories=%d timeout_s=%d",
            len(self.poi_categories),
            self.timeout,
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            self.area_buffered = self.area.buffer(0.008).simplify(0.005)

        assert self.area_buffered is not None
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

        log_progress(
            logger,
            self.verbose,
            "event=poi.query.build.done bbox=%s,%s,%s,%s duration_s=%.3f",
            lat_min,
            lon_min,
            lat_max,
            lon_max,
            time.perf_counter() - query_build_start,
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

        log_progress(
            logger,
            self.verbose,
            "event=poi.fetch.start endpoint=%s",
            request_header.removesuffix("?data="),
        )
        fetch_start = time.perf_counter()
        resp = requests.get(url=request_header + query_string, timeout=self.timeout)
        fetch_duration = time.perf_counter() - fetch_start
        log_progress(
            logger,
            self.verbose,
            "event=poi.fetch.response status_code=%d duration_s=%.3f",
            resp.status_code,
            fetch_duration,
        )

        if resp.status_code == 429:
            logger.warning(
                "event=poi.fetch.error status_code=429 duration_s=%.3f",
                fetch_duration,
            )
            raise RuntimeError(
                "429 Too Many Requests:\n"
                "You have sent multiple requests from the same IP and exceeded "
                "the fair use policy. Please wait a few minutes and try again."
            )
        elif resp.status_code == 504:
            logger.warning(
                "event=poi.fetch.error status_code=504 duration_s=%.3f",
                fetch_duration,
            )
            raise RuntimeError(
                "504 Gateway Timeout:\n"
                "The server is under heavy load and cannot process the request. "
                "Please try again later."
            )
        elif resp.status_code != 200:
            logger.error(
                "event=poi.fetch.error status_code=%d body=%s duration_s=%.3f",
                resp.status_code,
                resp.text[:160].replace("\n", " "),
                fetch_duration,
            )
            raise RuntimeError("Bad Request!")
        else:
            resp = json.loads(resp.text)

        parse_start = time.perf_counter()
        log_progress(
            logger,
            self.verbose,
            "event=poi.parse.start element_count=%d",
            len(resp["elements"]),
        )

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

        log_progress(
            logger,
            self.verbose,
            "event=poi.parse.split.done nodes=%d ways=%d",
            len(lst_nodes),
            len(lst_ways),
            level=logging.DEBUG,
        )

        nodes_df = pd.DataFrame(lst_nodes)
        ways_df = pd.DataFrame(lst_ways)

        if len(nodes_df) > 0 and len(ways_df) > 0:
            log_progress(
                logger,
                self.verbose,
                "event=poi.parse.join mode=nodes_and_ways",
                level=logging.DEBUG,
            )

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
            log_progress(
                logger,
                self.verbose,
                "event=poi.parse.join mode=ways_only",
                level=logging.DEBUG,
            )

            ways_df = ways_df.drop(columns=["id", "bounds", "nodes"])
            poi_df = ways_df.fillna(False)

        elif len(nodes_df) > 0 and len(ways_df) == 0:
            log_progress(
                logger,
                self.verbose,
                "event=poi.parse.join mode=nodes_only",
                level=logging.DEBUG,
            )

            nodes_df["geometry"] = nodes_df[["lon", "lat"]].apply(
                lambda p: [{"lat": p["lat"], "lon": p["lon"]}], axis=1
            )
            nodes_df = nodes_df.rename(
                columns={"lat": "center_latitude", "lon": "center_longitude"}
            )
            nodes_df = nodes_df.drop(columns=["id"])
            poi_df = nodes_df.fillna(False)
        else:
            logger.warning(
                "event=poi.parse.empty_response poi_categories=%d",
                len(self.poi_categories),
            )
            raise ValueError(
                "No POI data found for the specified poi_categories and area."
            )

        log_progress(
            logger,
            self.verbose,
            "event=poi.parse.done poi_count=%d duration_s=%.3f",
            len(poi_df),
            time.perf_counter() - parse_start,
        )

        for poi_category in self.poi_categories:
            if poi_category not in poi_df.columns:
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

        log_progress(
            logger,
            self.verbose,
            "event=poi.filter.done poi_count=%d",
            len(poi_df),
        )

        if len(poi_df) == 0:
            logger.warning("event=poi.filter.empty")
            raise ValueError("No POI data found within the study area.")

        return poi_df
