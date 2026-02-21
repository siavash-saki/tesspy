"""
Square and adaptive-square tessellation functions.
"""

import geopandas as gpd
import mercantile
import pandas as pd
from shapely.geometry import box


def get_squares_polyfill(gdf: gpd.GeoDataFrame, zoom_level: int) -> gpd.GeoDataFrame:
    """
    Square tessellation based on the quadKeys concept.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the area polygon(s)
    zoom_level : int
        Resolution, which controls the square sizes

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the squares
    """
    geom_name = gdf.geometry.name
    temp_dfs = []

    for _, rows in gdf.iterrows():
        gdf_geometry = rows[geom_name]
        bbox = gdf_geometry.bounds
        tiles = mercantile.tiles(bbox[0], bbox[1], bbox[2], bbox[3], zoom_level)
        temp_rows = []
        for tile in tiles:
            temp_row = rows.copy()
            square = box(*mercantile.bounds(tile))
            if square.intersects(gdf_geometry):
                temp_row[geom_name] = square
                temp_row["quadkey"] = mercantile.quadkey(tile)

                child_ids = mercantile.children(tile)
                temp_row["children_id"] = [
                    mercantile.quadkey(c_tile) for c_tile in child_ids
                ]

                temp_rows.append(temp_row)
        temp_dfs.append(pd.DataFrame(temp_rows))

    df = pd.concat(temp_dfs)
    df = df.reset_index(drop=True)

    gdf = gpd.GeoDataFrame(df, geometry=geom_name, crs="epsg:4326")

    return gdf


def get_adaptive_squares(
    input_gdf: gpd.GeoDataFrame, threshold: int
) -> gpd.GeoDataFrame:
    """
    Adaptive tessellation.

    Subdivides all squares where the POI count threshold is exceeded.

    Parameters
    ----------
    input_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles (polygons) with a 'count' column
    threshold : int
        Threshold, which controls the division of squares

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the updated squares
    """
    gdf = input_gdf.copy()
    gdf_exceeded = gdf[gdf["count"] >= threshold]

    for idx, row in gdf_exceeded.iterrows():
        children = gdf_exceeded.loc[[idx]]["children_id"].values[0]
        gdf.drop([idx], inplace=True)

        for child in children:
            new_row = row.copy()
            child_tile = mercantile.quadkey_to_tile(child)

            new_row["quadkey"] = child
            new_row["geometry"] = box(*mercantile.bounds(child_tile))
            grand_children = mercantile.children(child_tile)
            new_row["children_id"] = [
                mercantile.quadkey(c_tile) for c_tile in grand_children
            ]

            tmp_df = pd.DataFrame(new_row).transpose()
            tmp_gdf = gpd.GeoDataFrame(tmp_df, geometry="geometry", crs="epsg:4326")

            gdf = pd.concat([gdf, tmp_gdf], axis=0)

    gdf.index = gdf.reset_index(drop=True).index
    return gdf


def count_poi(df: gpd.GeoDataFrame, points: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Counts the number of POI in each tile.

    Parameters
    ----------
    df : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles (polygons)
    points : geopandas.GeoDataFrame
        GeoDataFrame containing the POI points

    Returns
    --------
    final_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles with an added 'count' column
    """
    pointsInPolygon = gpd.sjoin(df, points, how="left", predicate="contains")
    pointsInPolygon["count"] = 1
    pointsInPolygon.reset_index(inplace=True)

    tmp_a = pointsInPolygon.groupby(by="quadkey").count()
    tmp_a = tmp_a["count"].reset_index()
    tmp_a = tmp_a.sort_values(by="quadkey", ascending=True)

    tmp_b = df.reset_index().sort_values(by="quadkey", ascending=True)

    final_df = pd.merge(tmp_a, tmp_b, on="quadkey")
    final_gdf = final_df[["quadkey", "count", "geometry", "children_id"]]
    final_gdf = gpd.GeoDataFrame(final_gdf, geometry="geometry")

    return final_gdf
