"""
City block tessellation functions using road network data.
"""

import warnings

import geopandas as gpd
import numpy as np
from scipy.sparse import lil_matrix
from shapely import make_valid as shapely_make_valid
from shapely.geometry import LineString
from shapely.ops import polygonize, unary_union
from sklearn.cluster import AgglomerativeClustering

# ---------------------------------------------------------------------------
# New primary API
# ---------------------------------------------------------------------------


def create_city_blocks(
    road_network: gpd.GeoDataFrame,
    area: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Create city block polygons from a road network and study area boundary.

    The study area boundary is included as lines before polygonizing, so
    ``polygonize`` naturally produces closed blocks at the perimeter without
    requiring a separate gap-filling step.  All lines are noded via
    ``unary_union`` to split them at every intersection point.

    Parameters
    ----------
    road_network : geopandas.GeoDataFrame
        GeoDataFrame containing road LineString geometries.
    area : geopandas.GeoDataFrame
        GeoDataFrame with a single Polygon or MultiPolygon defining
        the study area boundary.  Must have a defined CRS.

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame of city block Polygons (CRS EPSG:4326).  Only blocks
        whose representative point falls within the study area are included.
    """
    area_polygon = area.geometry.iloc[0]

    # Extract boundary as individual LineStrings
    boundary = area_polygon.boundary
    if boundary.geom_type == "MultiLineString":
        boundary_lines = list(boundary.geoms)
    else:
        boundary_lines = [boundary]

    # Combine road geometries with boundary lines, then node everything
    # via unary_union (splits lines at every intersection point).
    all_lines = unary_union([*road_network.geometry, *boundary_lines])

    # Polygonize the fully-noded line network
    block_faces = list(polygonize(all_lines))

    if not block_faces:
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")

    blocks = gpd.GeoDataFrame(geometry=block_faces, crs="EPSG:4326")
    blocks["geometry"] = shapely_make_valid(blocks["geometry"].values)

    # Keep only blocks whose representative point is inside the study area.
    # representative_point() is guaranteed to lie inside the polygon.
    mask = blocks.geometry.representative_point().within(area_polygon)
    return blocks[mask].reset_index(drop=True)


def merge_city_blocks(
    blocks: gpd.GeoDataFrame,
    n_polygons: int,
) -> gpd.GeoDataFrame:
    """
    Merge city block polygons into *n_polygons* groups using
    adjacency-constrained agglomerative clustering.

    Only spatially adjacent blocks (those sharing a boundary segment) can
    be merged together, preventing disconnected MultiPolygon results.

    Parameters
    ----------
    blocks : geopandas.GeoDataFrame
        GeoDataFrame of city block Polygons.
    n_polygons : int
        Target number of merged blocks.

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame of merged block Polygons (CRS EPSG:4326).
    """
    n = len(blocks)

    # Build spatial adjacency matrix using the spatial index.
    adjacency = lil_matrix((n, n), dtype=np.int8)
    tree = blocks.sindex

    for i in range(n):
        candidates = tree.query(blocks.geometry.iloc[i], predicate="touches")
        for j in candidates:
            if j != i:
                adjacency[i, j] = 1
                adjacency[j, i] = 1

    connectivity = adjacency.tocsr()

    centroids = np.column_stack(
        [blocks.geometry.centroid.x, blocks.geometry.centroid.y]
    )

    model = AgglomerativeClustering(
        n_clusters=n_polygons,
        metric="euclidean",
        connectivity=connectivity,
    )
    model.fit(centroids)

    blocks = blocks.copy()
    blocks["_cluster"] = model.labels_
    blocks["geometry"] = shapely_make_valid(blocks["geometry"].values)

    merged = blocks.groupby("_cluster")["geometry"].agg(unary_union).tolist()
    result = gpd.GeoDataFrame({"geometry": merged}, crs=blocks.crs)

    # Explode any remaining MultiPolygons (rare with connectivity constraint)
    if (result.geom_type == "MultiPolygon").any():
        result = result.explode(index_parts=False).reset_index(drop=True)

    return result


# ---------------------------------------------------------------------------
# Deprecated helpers (kept for backward compatibility, will be removed v0.3.0)
# ---------------------------------------------------------------------------


def split_linestring(df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Split LineStrings with more than 2 points into 2-point segments.
    Each resulting segment retains the osmid of the original road.

    .. deprecated::
        Use :func:`create_city_blocks` instead, which handles noding
        automatically.  Will be removed in v0.3.0.

    Parameters
    ----------
    df : geopandas.GeoDataFrame
        GeoDataFrame containing shapely LineStrings

    Returns
    --------
    dataset : geopandas.GeoDataFrame
        GeoDataFrame with 2-point LineStrings
    """
    warnings.warn(
        "split_linestring is deprecated and will be removed in v0.3.0. "
        "Use create_city_blocks() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    linestrings = []
    osmid = []

    for geom, oid in zip(df.geometry, df["osmid"], strict=True):
        coords = geom.coords
        if len(coords) == 2:
            linestrings.append(geom)
            osmid.append(oid)
        else:
            for i in range(len(coords) - 1):
                linestrings.append(LineString(coords[i : i + 2]))
                osmid.append(oid)

    dataset = gpd.GeoDataFrame({"osmid": osmid, "geometry": linestrings})
    return dataset


def explode(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Explode MultiPolygon geometries into individual Polygon rows.

    .. deprecated::
        Use ``gdf.explode(index_parts=True)`` directly.
        Will be removed in v0.3.0.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame that may contain MultiPolygon geometries

    Returns
    --------
    gdf_out : geopandas.GeoDataFrame
        GeoDataFrame with only single Polygon geometries
    """
    warnings.warn(
        "explode is deprecated and will be removed in v0.3.0. "
        "Use gdf.explode(index_parts=True) directly.",
        DeprecationWarning,
        stacklevel=2,
    )
    return gdf.explode(index_parts=True)


def create_blocks(road_network: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Use shapely polygonize to create city block polygons from road LineStrings.

    .. deprecated::
        Use :func:`create_city_blocks` instead, which includes boundary
        handling and proper noding.  Will be removed in v0.3.0.

    Parameters
    ----------
    road_network : geopandas.GeoDataFrame
        GeoDataFrame containing street segment geometries

    Returns
    --------
    blocks : geopandas.GeoDataFrame
        GeoDataFrame of block polygons formed by the road network
    """
    warnings.warn(
        "create_blocks is deprecated and will be removed in v0.3.0. "
        "Use create_city_blocks() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    block_faces = list(polygonize(road_network["geometry"]))
    blocks = gpd.GeoDataFrame(geometry=block_faces).set_crs("EPSG:4326")
    return blocks


def get_rest_polygon(
    blocks: gpd.GeoDataFrame, area: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Create "rest polygons" to fill gaps not covered by road-based blocks.

    .. deprecated::
        Use :func:`create_city_blocks` instead, which includes the study
        area boundary in the polygonization so no gap-filling is needed.
        Will be removed in v0.3.0.

    Parameters
    ----------
    blocks : geopandas.GeoDataFrame
        GeoDataFrame containing city block polygons
    area : geopandas.GeoDataFrame
        GeoDataFrame containing the study area boundary polygon

    Returns
    --------
    rest_polygons : geopandas.GeoDataFrame
        GeoDataFrame containing the gap-filling polygons
    """
    warnings.warn(
        "get_rest_polygon is deprecated and will be removed in v0.3.0. "
        "Use create_city_blocks() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    blocks = blocks.copy()
    blocks["geometry"] = shapely_make_valid(blocks["geometry"].values)

    merged = unary_union(blocks["geometry"])

    rest = area.difference(gpd.GeoSeries([merged], crs="EPSG:4326"))
    rest = gpd.GeoDataFrame(rest, columns=["geometry"]).set_geometry("geometry")

    rest_polygons = rest.explode(index_parts=True).reset_index()
    rest_polygons = rest_polygons.drop(columns=["level_0"]).rename(
        columns={"level_1": "osm_id"}
    )

    return rest_polygons
