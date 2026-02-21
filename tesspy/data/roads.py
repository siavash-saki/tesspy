"""
Road network data retrieval via osmnx.
"""


import geopandas as gpd
import osmnx as ox

from tesspy._constants import OSM_HIGHWAY_TYPES


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
        If True, print progress information.
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

        if self.verbose:
            print(f"Selected highway type(s): {custom_filter}")

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
        if self.verbose:
            print("Collecting road network data...")
        graph = ox.graph_from_polygon(
            self.area.boundary.convex_hull.values[0], custom_filter=cf
        )
        graph_projected = ox.project_graph(graph, to_crs="epsg:4326")
        graph_undirected = graph_projected.to_undirected()
        graph_edges_as_gdf = ox.graph_to_gdfs(
            graph_undirected, nodes=False, edges=True
        )

        if self.verbose:
            print(f"Collected {len(graph_edges_as_gdf)} street segments.")

        return graph_edges_as_gdf
