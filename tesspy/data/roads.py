"""
Road network data retrieval via osmnx.
"""

import logging
import time

import geopandas as gpd
import osmnx as ox

from tesspy._constants import OSM_HIGHWAY_TYPES
from tesspy._logging import log_progress

logger = logging.getLogger(__name__)


class RoadData:
    """
    Retrieve and process OSM road network data using osmnx.

    Parameters
    ----------
    area : geopandas.GeoDataFrame
        GeoDataFrame with a single Polygon or MultiPolygon and a defined CRS.
    detail_deg : int or None
        Number of top highway types to include. None means all 19 types.
    split_roads : bool
        If True, LineStrings are split so each has exactly 2 points.
    verbose : bool
        If True, log progress information via the ``tesspy`` logger.
    """

    def __init__(
        self,
        area: gpd.GeoDataFrame,
        detail_deg: int | None = None,
        split_roads: bool = True,
        verbose: bool = False,
    ) -> None:
        self.detail_deg = detail_deg
        self.area = area
        self.verbose = verbose
        self.split_roads = split_roads

    @staticmethod
    def osm_highway_types() -> list[str]:
        """
        Return the list of OSM highway type categories.
        See https://wiki.openstreetmap.org/wiki/Key:highway

        Returns
        --------
        list of str
        """
        return OSM_HIGHWAY_TYPES

    def create_custom_filter(self) -> str:
        """
        Build the osmnx custom filter string for the selected highway types.

        Returns
        --------
        custom_filter : str
            Filter string in osmnx format, e.g. "['highway'~'motorway|trunk|...']"
        """
        if self.detail_deg is None:
            highwaytypes = self.osm_highway_types()
        elif isinstance(self.detail_deg, int):
            highwaytypes = self.osm_highway_types()[: self.detail_deg]
        else:
            raise ValueError("detail_deg must be None or an int")

        query = "|".join(highwaytypes[:-1]) + f"|{highwaytypes[-1]}"
        custom_filter = f"['highway'~'{query}']"

        log_progress(
            logger,
            self.verbose,
            "event=roads.filter.selected detail_deg=%s filter=%s",
            self.detail_deg,
            custom_filter,
        )

        return custom_filter

    def get_road_network(self) -> gpd.GeoDataFrame:
        """
        Download the road network for the study area and return it as a GeoDataFrame.

        Returns
        --------
        graph_edges_as_gdf : geopandas.GeoDataFrame
            GeoDataFrame containing road network edges
        """
        cf = self.create_custom_filter()
        log_progress(logger, self.verbose, "event=roads.fetch.start")
        fetch_start = time.perf_counter()
        graph = ox.graph_from_polygon(
            self.area.boundary.convex_hull.iloc[0], custom_filter=cf
        )
        graph_projected = ox.project_graph(graph, to_crs="epsg:4326")
        graph_undirected = graph_projected.to_undirected()
        graph_edges_as_gdf = ox.graph_to_gdfs(graph_undirected, nodes=False, edges=True)

        log_progress(
            logger,
            self.verbose,
            "event=roads.fetch.done segments=%d duration_s=%.3f",
            len(graph_edges_as_gdf),
            time.perf_counter() - fetch_start,
        )

        return graph_edges_as_gdf
