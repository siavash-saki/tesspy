"""
City block tessellation functions using road network data.
"""

import geopandas as gpd
from shapely.geometry import LineString, Point
from shapely.ops import polygonize, unary_union
from shapely.validation import make_valid


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

    for _, row in df.iterrows():
        if len(row["geometry"].coords) == 2:
            linestrings.append(row["geometry"])
            osmid.append(row["osmid"])
        else:
            for i in range(0, len(row["geometry"].coords) - 1):
                p1 = Point(row["geometry"].coords[i][0], row["geometry"].coords[i][1])
                p2 = Point(
                    row["geometry"].coords[i + 1][0], row["geometry"].coords[i + 1][1]
                )
                linestrings.append(LineString([p1, p2]))
                osmid.append(row["osmid"])

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
    gs = gdf.explode(index_parts=True)
    gdf2 = gs.reset_index().rename(columns={0: "geometry"})
    gdf_out = gdf2.merge(
        gdf.drop("geometry", axis=1),
        left_on="level_0",
        right_index=True,
    )
    gdf_out = gdf_out.set_index(["level_0", "level_1"]).set_geometry("geometry")
    gdf_out.crs = gdf.crs
    return gdf_out


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
        blocks["geometry"] = blocks["geometry"].apply(lambda x: make_valid(x))

        merged_polygons = gpd.GeoSeries(unary_union(blocks["geometry"].values))
        merged_polygons.set_crs("EPSG:4326", allow_override=True, inplace=True)

        rest = area.difference(merged_polygons)
        rest = gpd.GeoDataFrame(rest)
        rest = rest.rename(columns={0: "geometry"}).set_geometry("geometry")

        rest_polygons = explode(rest)
        rest_polygons.reset_index(inplace=True)
        rest_polygons.drop(columns=["level_0"], inplace=True)
        rest_polygons.rename(columns={"level_1": "osm_id"}, inplace=True)

        return rest_polygons

    else:
        raise ValueError("City blocks and the area both require a geometry attribute.")
