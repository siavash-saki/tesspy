"""
POI (Point of Interest) data retrieval via osmnx (Overpass API).
"""

import logging
import time

import osmnx as ox
import pandas as pd

from tesspy._constants import OSM_PRIMARY_FEATURES
from tesspy._logging import log_progress

logger = logging.getLogger(__name__)


class POIdata:
    """
    Query OSM for Points of Interest within a study area using osmnx.

    osmnx provides automatic retry on 429/504 errors, proactive rate-limit
    checking, file-based response caching, and polygon-based spatial
    filtering.

    Parameters
    ----------
    area : geopandas.GeoDataFrame
        GeoDataFrame with a single Polygon or MultiPolygon and a defined CRS.
    poi_categories : list of str
        OSM primary map feature categories to query.
    timeout : int
        Request timeout in seconds for the Overpass query.
    verbose : bool
        If True, log progress information via the ``tesspy`` logger.
    """

    def __init__(
        self,
        area,
        poi_categories: list[str],
        timeout: int,
        verbose: bool,
    ) -> None:
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

    def get_poi_data(self) -> pd.DataFrame:
        """
        Query OSM POI data and return it as a DataFrame.

        Returns
        --------
        poi_df : pandas.DataFrame
            DataFrame with center coordinates and boolean columns for each
            queried POI category.
        """
        # Validate categories
        for poi_category in self.poi_categories:
            if poi_category not in self.osm_primary_features():
                raise ValueError(
                    f"{poi_category} is not a valid POI primary category. "
                    f"See a list of OSM primary features with "
                    f"Tessellation.osm_primary_features()"
                )

        tags: dict[str, bool | str | list[str]] = dict.fromkeys(
            self.poi_categories, True
        )
        polygon = self.area.geometry.iloc[0]

        log_progress(
            logger,
            self.verbose,
            "event=poi.fetch.start poi_categories=%d timeout_s=%d",
            len(self.poi_categories),
            self.timeout,
        )

        # Configure osmnx timeout for this query
        ox.settings.requests_timeout = self.timeout

        fetch_start = time.perf_counter()
        gdf = ox.features_from_polygon(polygon, tags=tags)
        fetch_duration = time.perf_counter() - fetch_start

        log_progress(
            logger,
            self.verbose,
            "event=poi.fetch.done features=%d duration_s=%.3f",
            len(gdf),
            fetch_duration,
        )

        if len(gdf) == 0:
            logger.warning(
                "event=poi.fetch.empty poi_categories=%d",
                len(self.poi_categories),
            )
            raise ValueError(
                "No POI data found for the specified poi_categories and area."
            )

        # Extract centroid coordinates
        centroids = gdf.geometry.centroid
        result = pd.DataFrame(
            {
                "center_longitude": centroids.x.values,
                "center_latitude": centroids.y.values,
            }
        )

        # Add boolean category columns
        for cat in self.poi_categories:
            if cat in gdf.columns:
                result[cat] = gdf[cat].notna().values
            else:
                result[cat] = False

        result = result.reset_index(drop=True)

        log_progress(
            logger,
            self.verbose,
            "event=poi.parse.done poi_count=%d",
            len(result),
        )

        return result
