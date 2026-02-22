"""
City block tessellation functions using road network data.
"""

import geopandas as gpd
from shapely import make_valid as shapely_make_valid
from shapely.geometry import LineString
from shapely.ops import polygonize, unary_union


def split_linestring(df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Split LineStrings with more than 2 points into 2-point segments.
    Each resulting segment retains the osmid of the original road.

    Parameters
    ----------
    df : geopandas.GeoDataFrame
        GeoDataFrame containing shapely LineStrings

    Returns
    --------
    dataset : geopandas.GeoDataFrame
        GeoDataFrame with 2-point LineStrings
    """
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

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame that may contain MultiPolygon geometries

    Returns
    --------
    gdf_out : geopandas.GeoDataFrame
        GeoDataFrame with only single Polygon geometries
    """
    return gdf.explode(index_parts=True)


def create_blocks(road_network: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Use shapely polygonize to create city block polygons from road LineStrings.

    Parameters
    ----------
    road_network : geopandas.GeoDataFrame
        GeoDataFrame containing street segment geometries

    Returns
    --------
    blocks : geopandas.GeoDataFrame
        GeoDataFrame of block polygons formed by the road network
    """
    if hasattr(road_network, "geometry"):
        block_faces = list(polygonize(road_network["geometry"]))
        blocks = gpd.GeoDataFrame(geometry=block_faces).set_crs("EPSG:4326")
        return blocks
    else:
        raise AttributeError("Road network data must have a geometry attribute.")


def get_rest_polygon(
    blocks: gpd.GeoDataFrame, area: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Create "rest polygons" to fill gaps not covered by road-based blocks.

    Dead-ends and boundary regions that cannot form closed blocks are filled
    by subtracting the union of all blocks from the study area boundary.

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
    if hasattr(blocks, "geometry") and hasattr(area, "geometry"):
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

    else:
        raise ValueError("City blocks and the area both require a geometry attribute.")
