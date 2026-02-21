"""
Core Tessellation class — the primary public interface of tesspy.
"""

import logging
import warnings
from typing import Literal

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi
from shapely.geometry import MultiPolygon, Point
from sklearn.cluster import KMeans

import hdbscan
from sklearn.cluster import AgglomerativeClustering
from shapely.ops import unary_union

from tesspy._constants import (
    DEFAULT_POI_CATEGORIES,
    OSM_HIGHWAY_TYPES,
    OSM_PRIMARY_FEATURES,
)
from tesspy._validators import _check_input_geodataframe, _check_valid_geometry_gdf
from tesspy.data._geo import count_poi_per_tile, get_city_polygon
from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData
from tesspy.methods.city_blocks import (
    create_blocks,
    explode,
    get_rest_polygon,
    split_linestring,
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
        return [
            cat for cat in poi_list if cat not in self.poi_dataframe.columns
        ]

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
        poi_categories: list[str] | Literal["all"] = None,
        threshold: int | None = None,
        timeout: int = 60,
        verbose: bool = False,
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
        verbose : bool, default=False
            Log progress information via the ``tesspy`` logger

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with adaptive square tiles
        """
        if poi_categories is None:
            poi_categories = DEFAULT_POI_CATEGORIES.copy()

        if poi_categories == "all":
            poi_categories = self.osm_primary_features()

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
                logger.info(
                    "Threshold=%d  => set as the median POI count "
                    "per square at the initial level",
                    threshold,
                )

        i = start_resolution
        while max(aqk_count_df["count"].values) > threshold:
            i += 1
            if verbose:
                logger.info(
                    "Threshold exceeded. Subdividing to resolution %d...",
                    i,
                )

            df_tmp = get_adaptive_squares(aqk_count_df, threshold)
            df_tmp.drop(columns=["count"], inplace=True)
            df_tmp2 = count_poi(df_tmp, poi_data_aqk)
            aqk_count_df = df_tmp2

        final_aqk = gpd.sjoin(aqk_count_df, self.area_gdf)
        final_aqk = final_aqk.drop(columns=["osm_id", "children_id", "index_right"])

        return final_aqk

    def voronoi(
        self,
        cluster_algo: Literal["k-means", "hdbscan"] | None = "k-means",
        poi_categories: list[str] | Literal["all"] = None,
        timeout: int = 60,
        n_polygons: int = 100,
        min_cluster_size: int = 15,
        verbose: bool = False,
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
        verbose : bool, default=False
            Log progress information via the ``tesspy`` logger

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with Voronoi polygon tiles and a 'voronoi_id' column
        """
        if poi_categories is None:
            poi_categories = DEFAULT_POI_CATEGORIES.copy()

        if poi_categories == "all":
            poi_categories = self.osm_primary_features()

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

        tess_data = self.poi_dataframe[
            self.poi_dataframe[poi_categories].sum(axis=1) > 0
        ]
        data_locs = tess_data[["center_longitude", "center_latitude"]].values

        if cluster_algo == "k-means":
            if verbose:
                logger.info("K-Means Clustering...")
            clustering = KMeans(n_clusters=n_polygons).fit(data_locs)
            generators = [
                np.mean(data_locs[clustering.labels_ == label], axis=0)
                for label in range(n_polygons)
            ]

        elif cluster_algo == "hdbscan":
            if verbose:
                logger.info("HDBSCAN Clustering... This can take a while...")
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
                    "Too many generators for Voronoi diagram. "
                    "Please select a clustering algorithm."
                )
            else:
                generators = data_locs

        else:
            raise ValueError(
                "cluster_algo must be one of: 'k-means', 'hdbscan', or None"
            )

        if verbose:
            logger.info("Creating Voronoi polygons...")
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
        self,
        n_polygons: int | None = None,
        detail_deg: int | None = None,
        split_roads: bool = True,
        verbose: bool = False,
    ) -> gpd.GeoDataFrame:
        """
        Create city block tiles using OSM road network data.

        Parameters
        ----------
        n_polygons : int or None, default=None
            Target number of city blocks (approximate). Uses hierarchical
            clustering to merge small blocks. If None, all raw blocks are returned.
        detail_deg : int or None, default=None
            Number of top OSM highway types to include. None means all 19 types.
        split_roads : bool, default=True
            Split LineStrings so each has exactly 2 points (more robust
            polygonization, but slower).
        verbose : bool, default=False
            Log progress information via the ``tesspy`` logger

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with city block tiles and a 'cityblock_id' column
        """
        if detail_deg is None:
            highwaytypes = self.osm_highway_types()
        elif (
            isinstance(detail_deg, int)
            and detail_deg <= len(self.osm_highway_types())
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
                queried_area, detail_deg, split_roads, verbose
            )
            road_data = road_data_collect_object.get_road_network()
            self.road_network = road_data
            self.queried_highway_types = highwaytypes
        else:
            road_data = self.road_network

        if split_roads:
            if verbose:
                logger.info("Splitting LineStrings to 2-point segments...")
            road_data = split_linestring(road_data)

        if verbose:
            logger.info("Creating initial city blocks from road network...")

        blocks = create_blocks(road_data)

        polygons_in_area = gpd.sjoin(blocks, queried_area, how="inner")
        polygons_in_area.drop(columns=["index_right"], inplace=True)

        rest_polygons = get_rest_polygon(polygons_in_area, queried_area)

        city_blocks = pd.concat([polygons_in_area, rest_polygons])

        if not n_polygons:
            city_blocks = city_blocks[["geometry"]].reset_index(drop=True)
            city_blocks = _check_valid_geometry_gdf(city_blocks)
            city_blocks.reset_index(inplace=True)
            city_blocks.rename(columns={"index": "cityblock_id"}, inplace=True)
            city_blocks["cityblock_id"] = (
                "cityblockID" + city_blocks["cityblock_id"].astype(str)
            )
            return city_blocks

        if n_polygons > len(city_blocks):
            raise ValueError(
                f"Cannot extract more city blocks than the initial count. "
                f"Initial: {len(city_blocks)}, requested: {n_polygons}. "
                f"Choose a value less than {len(city_blocks)}."
            )

        if verbose:
            logger.info("Merging small city blocks with hierarchical clustering...")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            city_blocks["centroid"] = city_blocks.centroid

        coordinates = np.column_stack(
            [city_blocks["centroid"].x, city_blocks["centroid"].y]
        )
        # Note: AgglomerativeClustering requires O(n²) memory and O(n³) time.
        # Large datasets may exhaust RAM.
        model = AgglomerativeClustering(n_clusters=n_polygons, metric="euclidean")
        model.fit(coordinates)

        city_blocks["Cluster"] = model.labels_

        merged_polys = []
        for idx in city_blocks["Cluster"].unique():
            tmp = city_blocks[city_blocks["Cluster"] == idx]
            polygons = tmp["geometry"].to_numpy()
            merged_polygon = gpd.GeoSeries(unary_union(polygons))
            merged_polys.append(merged_polygon[0])

        merged_polys_df = gpd.GeoDataFrame({"geometry": merged_polys}, crs="EPSG:4326")
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
        final_city_blocks["cityblock_id"] = (
            "cityblockID" + final_city_blocks["cityblock_id"].astype(str)
        )

        return final_city_blocks

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
