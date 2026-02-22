"""
Core Tessellation class — the primary public interface of tesspy.
"""

import logging
import time
from typing import Literal

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi
from shapely.geometry import MultiPolygon
from sklearn.cluster import HDBSCAN, KMeans

from tesspy._constants import (
    DEFAULT_POI_CATEGORIES,
    OSM_HIGHWAY_TYPES,
    OSM_PRIMARY_FEATURES,
)
from tesspy._logging import log_progress
from tesspy._validators import _check_input_geodataframe, _check_valid_geometry_gdf
from tesspy.data._geo import count_poi_per_tile, get_city_polygon
from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData
from tesspy.methods.city_blocks import (
    create_city_blocks,
    merge_city_blocks,
)
from tesspy.methods.hexagons import get_h3_hexagons
from tesspy.methods.squares import count_poi, get_adaptive_squares, get_squares_polyfill
from tesspy.methods.voronoi import voronoi_polygons

logger = logging.getLogger(__name__)

# Backward-compatibility re-exports so that code doing
#   from tesspy.tessellation import get_city_polygon
#   from tesspy.tessellation import count_poi_per_tile
# continues to work.
__all__ = [
    "Tessellation",
    "get_city_polygon",
    "count_poi_per_tile",
    "_check_input_geodataframe",
    "_check_valid_geometry_gdf",
]


class Tessellation:
    """
    Create a Tessellation object for a geographic area, enabling multiple
    tessellation methods (squares, hexagons, adaptive squares, Voronoi,
    and city blocks).

    Parameters
    ----------
    area : geopandas.GeoDataFrame or str
        GeoDataFrame must have a single Polygon or MultiPolygon in its
        geometry column with a defined CRS.
        str must be a city name or address of a region.

    Examples
    --------
    >>> ffm = Tessellation('Frankfurt am Main')
    >>> squares_gdf = ffm.squares(resolution=12)
    """

    def __init__(self, area: gpd.GeoDataFrame | str) -> None:
        if isinstance(area, gpd.GeoDataFrame):
            self.area_gdf = _check_input_geodataframe(area)
        elif isinstance(area, str):
            self.area_gdf = get_city_polygon(area)
        else:
            raise TypeError("area must be a GeoDataFrame or a string (city name)")

        self.poi_dataframe: pd.DataFrame = pd.DataFrame()
        self.available_poi_categories: list[str] = []
        self.road_network: pd.DataFrame = pd.DataFrame()
        self.queried_highway_types: list[str] = []

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _get_missing_poi_categories(self, poi_list: list[str]) -> list[str]:
        """Return categories in poi_list that are not yet in self.poi_dataframe."""
        return [cat for cat in poi_list if cat not in self.poi_dataframe.columns]

    # ------------------------------------------------------------------
    # Tessellation methods
    # ------------------------------------------------------------------

    def squares(self, resolution: int) -> gpd.GeoDataFrame:
        """
        Generate a regular square grid over the area.

        Parameters
        ----------
        resolution : int
            Zoom level controlling square size (higher = smaller squares)

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame containing the square tiles
        """
        df_qk_squares = get_squares_polyfill(self.area_gdf, resolution)
        df_qk_squares = df_qk_squares.drop(columns=["osm_id", "children_id"])
        return df_qk_squares

    def hexagons(self, resolution: int) -> gpd.GeoDataFrame:
        """
        Generate a regular hexagon grid over the area (Uber H3).

        Parameters
        ----------
        resolution : int
            H3 resolution controlling hexagon size (0–15)

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame containing the hexagon tiles with a 'hex_id' column
        """
        df_h3_hexagons = get_h3_hexagons(self.area_gdf, resolution)
        df_h3_hexagons = df_h3_hexagons.reset_index().rename(
            columns={"index": "hex_id"}
        )
        return df_h3_hexagons

    def adaptive_squares(
        self,
        start_resolution: int,
        poi_categories: list[str] | Literal["all"] | None = None,
        threshold: int | None = None,
        timeout: int = 60,
        verbose: bool = True,
    ) -> gpd.GeoDataFrame:
        """
        Generate adaptive squares based on POI density.

        Squares are created at start_resolution and recursively subdivided
        into four smaller squares whenever the POI count exceeds the threshold.

        Parameters
        ----------
        start_resolution : int
            Initial zoom level for the square grid
        poi_categories : list of str or 'all', default=["amenity", "building"]
            OSM primary map feature categories used to guide subdivision
        threshold : int or None, default=None
            POI count threshold for subdivision. If None, uses the median
            count across all initial squares.
        timeout : int, default=60
            Overpass API timeout in seconds
        verbose : bool, default=True
            Log progress information via the ``tesspy`` logger

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with adaptive square tiles
        """
        method_start = time.perf_counter()
        if poi_categories is None:
            poi_categories = DEFAULT_POI_CATEGORIES.copy()

        if poi_categories == "all":
            poi_categories = self.osm_primary_features()

        log_progress(
            logger,
            verbose,
            "event=adaptive_squares.start start_resolution=%d poi_categories=%d "
            "threshold=%s timeout_s=%d",
            start_resolution,
            len(poi_categories),
            threshold,
            timeout,
        )

        missing_poi_categories = self._get_missing_poi_categories(poi_categories)

        if len(missing_poi_categories) > 0:
            poi_data_obj = POIdata(
                self.area_gdf, missing_poi_categories, timeout, verbose
            )
            poi_data_new = poi_data_obj.get_poi_data()
            self.poi_dataframe = pd.concat([self.poi_dataframe, poi_data_new]).fillna(
                False
            )
            self.poi_dataframe = self.poi_dataframe.reset_index(drop=True)
        else:
            log_progress(
                logger,
                verbose,
                "event=poi.cache.hit categories=%d",
                len(poi_categories),
                level=logging.DEBUG,
            )

        tess_data = self.poi_dataframe[
            self.poi_dataframe[poi_categories].sum(axis=1) > 0
        ]
        points_geom = gpd.points_from_xy(
            tess_data["center_longitude"], tess_data["center_latitude"]
        )
        tess_data = gpd.GeoDataFrame(
            geometry=points_geom, data=tess_data[poi_categories], crs="EPSG:4326"
        )
        poi_data_aqk = tess_data.rename(columns={"points_geom": "geometry"})
        df_aqk = get_squares_polyfill(self.area_gdf, start_resolution)
        aqk_count_df = count_poi(df_aqk, poi_data_aqk)

        if not threshold:
            threshold = int(aqk_count_df["count"].median())
            log_progress(
                logger,
                verbose,
                "event=adaptive_squares.threshold.auto threshold=%d",
                threshold,
            )

        i = start_resolution
        while aqk_count_df["count"].max() > threshold:
            i += 1
            log_progress(
                logger,
                verbose,
                "event=adaptive_squares.subdivide resolution=%d",
                i,
            )

            df_tmp = get_adaptive_squares(aqk_count_df, threshold)
            df_tmp = df_tmp.drop(columns=["count"])
            aqk_count_df = count_poi(df_tmp, poi_data_aqk)

        final_aqk = gpd.sjoin(aqk_count_df, self.area_gdf)
        final_aqk = final_aqk.drop(columns=["osm_id", "children_id", "index_right"])

        log_progress(
            logger,
            verbose,
            "event=adaptive_squares.done tiles=%d duration_s=%.3f",
            len(final_aqk),
            time.perf_counter() - method_start,
        )

        return final_aqk

    def voronoi(
        self,
        cluster_algo: Literal["k-means", "hdbscan"] | None = "k-means",
        poi_categories: list[str] | Literal["all"] | None = None,
        timeout: int = 60,
        n_polygons: int = 100,
        min_cluster_size: int = 15,
        verbose: bool = True,
    ) -> gpd.GeoDataFrame:
        """
        Generate Voronoi polygon tessellation driven by POI density.

        Parameters
        ----------
        cluster_algo : {'k-means', 'hdbscan', None}, default='k-means'
            Clustering algorithm used to derive Voronoi generators.
            If None, POI locations are used directly (max 5000).
        poi_categories : list of str or 'all', default=["amenity", "building"]
            OSM primary map feature categories used as input data
        timeout : int, default=60
            Overpass API timeout in seconds
        n_polygons : int, default=100
            Target number of polygons (k-means only)
        min_cluster_size : int, default=15
            Minimum cluster size (hdbscan only)
        verbose : bool, default=True
            Log progress information via the ``tesspy`` logger

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with Voronoi polygon tiles and a 'voronoi_id' column
        """
        method_start = time.perf_counter()
        if poi_categories is None:
            poi_categories = DEFAULT_POI_CATEGORIES.copy()

        if poi_categories == "all":
            poi_categories = self.osm_primary_features()

        log_progress(
            logger,
            verbose,
            "event=voronoi.start cluster_algo=%s poi_categories=%d n_polygons=%d "
            "min_cluster_size=%d timeout_s=%d",
            cluster_algo,
            len(poi_categories),
            n_polygons,
            min_cluster_size,
            timeout,
        )

        missing_poi_categories = self._get_missing_poi_categories(poi_categories)

        if isinstance(self.area_gdf.geometry.iloc[0], MultiPolygon):
            queried_area = self.area_gdf.convex_hull
        else:
            queried_area = self.area_gdf

        if len(missing_poi_categories) > 0:
            poi_data_obj = POIdata(
                queried_area, missing_poi_categories, timeout, verbose
            )
            poi_data_new = poi_data_obj.get_poi_data()
            self.poi_dataframe = pd.concat([self.poi_dataframe, poi_data_new]).fillna(
                False
            )
            self.poi_dataframe = self.poi_dataframe.reset_index(drop=True)
        else:
            log_progress(
                logger,
                verbose,
                "event=poi.cache.hit categories=%d",
                len(poi_categories),
                level=logging.DEBUG,
            )

        tess_data = self.poi_dataframe[
            self.poi_dataframe[poi_categories].sum(axis=1) > 0
        ]
        data_locs = tess_data[["center_longitude", "center_latitude"]].to_numpy()

        if cluster_algo == "k-means":
            log_progress(logger, verbose, "event=voronoi.cluster.start algo=kmeans")
            clustering = KMeans(n_clusters=n_polygons).fit(data_locs)
            generators = np.empty((n_polygons, 2))
            for label in range(n_polygons):
                generators[label] = data_locs[clustering.labels_ == label].mean(axis=0)

        elif cluster_algo == "hdbscan":
            log_progress(logger, verbose, "event=voronoi.cluster.start algo=hdbscan")
            clustering = HDBSCAN(
                min_cluster_size=min_cluster_size,
            ).fit(data_locs)
            n_clusters = clustering.labels_.max() + 1
            generators = np.empty((n_clusters, 2))
            for label in range(n_clusters):
                generators[label] = data_locs[clustering.labels_ == label].mean(axis=0)

        elif cluster_algo is None:
            if len(tess_data) > 5000:
                raise ValueError(
                    "Too many generators for Voronoi diagram. "
                    "Please select a clustering algorithm."
                )
            else:
                generators = data_locs

        else:
            raise ValueError(
                "cluster_algo must be one of: 'k-means', 'hdbscan', or None"
            )

        log_progress(
            logger,
            verbose,
            "event=voronoi.polygons.create generators=%d",
            len(generators),
        )
        voronoi_dia = Voronoi(generators)
        voronoi_poly = gpd.GeoDataFrame(
            geometry=list(voronoi_polygons(voronoi_dia, 0.1)), crs="EPSG:4326"
        )
        voronoi_poly = gpd.sjoin(voronoi_poly, self.area_gdf)
        vor_polygons = voronoi_poly.intersection(self.area_gdf.geometry.iloc[0])
        df_voronoi = gpd.GeoDataFrame(geometry=vor_polygons)

        df_voronoi = _check_valid_geometry_gdf(df_voronoi)

        df_voronoi = df_voronoi.reset_index()
        df_voronoi = df_voronoi.rename(columns={"index": "voronoi_id"})
        df_voronoi["voronoi_id"] = "voronoiID" + df_voronoi["voronoi_id"].astype(str)

        log_progress(
            logger,
            verbose,
            "event=voronoi.done polygons=%d duration_s=%.3f",
            len(df_voronoi),
            time.perf_counter() - method_start,
        )

        return df_voronoi

    def city_blocks(
        self,
        n_polygons: int | None = None,
        detail_deg: int | None = None,
        verbose: bool = True,
    ) -> gpd.GeoDataFrame:
        """
        Create city block tiles using OSM road network data.

        The study area boundary is included in the road line network before
        polygonization, producing gap-free blocks that fully tile the area.
        When *n_polygons* is set, adjacent blocks are merged using
        connectivity-constrained agglomerative clustering.

        Parameters
        ----------
        n_polygons : int or None, default=None
            Target number of city blocks (approximate).  Uses adjacency-
            constrained hierarchical clustering to merge small blocks.
            If None, all raw blocks are returned.
        detail_deg : int or None, default=None
            Number of top OSM highway types to include. None means all 19 types.
        verbose : bool, default=True
            Log progress information via the ``tesspy`` logger

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with city block tiles and a 'cityblock_id' column
        """
        method_start = time.perf_counter()
        log_progress(
            logger,
            verbose,
            "event=city_blocks.start n_polygons=%s detail_deg=%s",
            n_polygons,
            detail_deg,
        )

        if detail_deg is None:
            highwaytypes = self.osm_highway_types()
        elif isinstance(detail_deg, int) and detail_deg <= len(
            self.osm_highway_types()
        ):
            highwaytypes = self.osm_highway_types()[:detail_deg]
        else:
            raise ValueError("detail_deg must be None or a valid int")

        if isinstance(self.area_gdf.geometry.iloc[0], MultiPolygon):
            queried_area = self.area_gdf.convex_hull
        else:
            queried_area = self.area_gdf

        if self.queried_highway_types != highwaytypes:
            road_data_collect_object = RoadData(
                queried_area, detail_deg, verbose=verbose
            )
            road_data = road_data_collect_object.get_road_network()
            self.road_network = road_data
            self.queried_highway_types = highwaytypes
        else:
            road_data = self.road_network
            log_progress(
                logger,
                verbose,
                "event=roads.cache.hit highwaytypes=%d",
                len(highwaytypes),
                level=logging.DEBUG,
            )

        log_progress(logger, verbose, "event=city_blocks.blocks.create.start")

        city_blocks_gdf = create_city_blocks(road_data, self.area_gdf)

        if n_polygons is not None:
            if n_polygons > len(city_blocks_gdf):
                raise ValueError(
                    f"Cannot extract more city blocks than the initial count. "
                    f"Initial: {len(city_blocks_gdf)}, requested: {n_polygons}. "
                    f"Choose a value less than {len(city_blocks_gdf)}."
                )
            log_progress(logger, verbose, "event=city_blocks.merge.start")
            city_blocks_gdf = merge_city_blocks(city_blocks_gdf, n_polygons)

        city_blocks_gdf = _check_valid_geometry_gdf(city_blocks_gdf)
        city_blocks_gdf = city_blocks_gdf[["geometry"]].reset_index(drop=True)
        city_blocks_gdf = city_blocks_gdf.reset_index()
        city_blocks_gdf = city_blocks_gdf.rename(columns={"index": "cityblock_id"})
        city_blocks_gdf["cityblock_id"] = (
            "cityblockID" + city_blocks_gdf["cityblock_id"].astype(str)
        )

        log_progress(
            logger,
            verbose,
            "event=city_blocks.done polygons=%d duration_s=%.3f",
            len(city_blocks_gdf),
            time.perf_counter() - method_start,
        )

        return city_blocks_gdf

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get_polygon(self) -> gpd.GeoDataFrame:
        """Return the study area polygon as a GeoDataFrame."""
        return self.area_gdf

    def get_poi_data(self) -> pd.DataFrame:
        """Return the cached POI DataFrame (empty if no POI methods called yet)."""
        return self.poi_dataframe

    def get_road_network(self) -> pd.DataFrame:
        """Return the cached road network GeoDataFrame."""
        return self.road_network

    # ------------------------------------------------------------------
    # Static helpers
    # ------------------------------------------------------------------

    @staticmethod
    def osm_primary_features() -> list[str]:
        """
        Return the list of primary OSM map feature categories.
        See https://wiki.openstreetmap.org/wiki/Map_features
        """
        return OSM_PRIMARY_FEATURES

    @staticmethod
    def osm_highway_types() -> list[str]:
        """
        Return the list of OSM highway type categories.
        See https://wiki.openstreetmap.org/wiki/Key:highway
        """
        return OSM_HIGHWAY_TYPES
