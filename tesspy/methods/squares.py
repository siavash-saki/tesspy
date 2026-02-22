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
    records: list[dict] = []

    for _, rows in gdf.iterrows():
        gdf_geometry = rows[geom_name]
        row_data = {k: v for k, v in rows.items() if k != geom_name}
        bbox = gdf_geometry.bounds
        tiles = mercantile.tiles(bbox[0], bbox[1], bbox[2], bbox[3], zoom_level)
        for tile in tiles:
            square = box(*mercantile.bounds(tile))
            if square.intersects(gdf_geometry):
                record = row_data.copy()
                record[geom_name] = square
                record["quadkey"] = mercantile.quadkey(tile)
                child_ids = mercantile.children(tile)
                record["children_id"] = [
                    mercantile.quadkey(c_tile) for c_tile in child_ids
                ]
                records.append(record)

    gdf = gpd.GeoDataFrame(records, geometry=geom_name, crs="epsg:4326")

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
    exceeded_mask = input_gdf["count"] >= threshold
    kept = input_gdf[~exceeded_mask]
    exceeded = input_gdf[exceeded_mask]

    new_rows: list[dict] = []
    for children in exceeded["children_id"]:
        for child in children:
            child_tile = mercantile.quadkey_to_tile(child)
            grand_children = mercantile.children(child_tile)
            new_rows.append(
                {
                    "quadkey": child,
                    "geometry": box(*mercantile.bounds(child_tile)),
                    "children_id": [
                        mercantile.quadkey(c_tile) for c_tile in grand_children
                    ],
                }
            )

    if new_rows:
        new_gdf = gpd.GeoDataFrame(new_rows, geometry="geometry", crs="epsg:4326")
        gdf = pd.concat([kept, new_gdf], ignore_index=True)
    else:
        gdf = kept.reset_index(drop=True)

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
    joined = gpd.sjoin(df, points, how="left", predicate="contains")
    counts = joined.groupby("quadkey").size().rename("count")

    final_gdf = df[["quadkey", "geometry", "children_id"]].merge(
        counts, on="quadkey", how="left"
    )
    final_gdf["count"] = final_gdf["count"].fillna(0).astype(int)
    final_gdf = gpd.GeoDataFrame(final_gdf, geometry="geometry")

    return final_gdf
