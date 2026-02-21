"""
Hexagon tessellation functions (Uber H3).
"""

import geopandas as gpd
import h3
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon


def get_h3_hexagons(gdf: gpd.GeoDataFrame, resolution: int) -> gpd.GeoDataFrame:
    """
    Hexagon tessellation based on the H3 implementation by Uber.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the area polygon(s)
    resolution : int
        Resolution, which controls the hexagon sizes

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the hexagons
    """
    if isinstance(gdf.geometry.iloc[0], Polygon):
        hexs = h3.polyfill(
            gdf.geometry[0].__geo_interface__, resolution, geo_json_conformant=True
        )
        polygonise = lambda hex_id: Polygon(
            h3.h3_to_geo_boundary(hex_id, geo_json=True)
        )
        all_polys = gpd.GeoSeries(
            list(map(polygonise, hexs)), index=hexs, crs="EPSG:4326"
        )

        gdf = gpd.GeoDataFrame(geometry=all_polys, crs="EPSG:4326")
        return gdf

    elif isinstance(gdf.geometry.iloc[0], MultiPolygon):
        parts_lst = []
        for idx, row in gdf.explode(index_parts=True).loc[0].iterrows():
            hexs = h3.polyfill(
                row.geometry.__geo_interface__, resolution, geo_json_conformant=True
            )

            polygonise = lambda hex_id: Polygon(
                h3.h3_to_geo_boundary(hex_id, geo_json=True)
            )
            all_polys = gpd.GeoSeries(
                list(map(polygonise, hexs)), index=hexs, crs="EPSG:4326"
            )

            part_gdf = gpd.GeoDataFrame(geometry=all_polys)
            parts_lst.append(part_gdf)

        return pd.concat(parts_lst)
