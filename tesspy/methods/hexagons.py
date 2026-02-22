"""
Hexagon tessellation functions (Uber H3).
"""

import geopandas as gpd
import h3
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon


def _hex_to_polygon(hex_id: str) -> Polygon:
    """Convert an H3 hex ID to a shapely Polygon."""
    # cell_to_boundary returns (lat, lng) tuples; shapely needs (lng, lat)
    coords = h3.cell_to_boundary(hex_id)
    return Polygon([(lng, lat) for lat, lng in coords])


def _geojson_to_h3poly(geojson: dict) -> h3.LatLngPoly:
    """Convert a GeoJSON Polygon dict to an H3Poly.

    GeoJSON uses (lng, lat) order; H3 uses (lat, lng).
    """
    coords = geojson["coordinates"]
    outer = [(lat, lng) for lng, lat in coords[0]]
    holes = [[(lat, lng) for lng, lat in hole] for hole in coords[1:]]
    return h3.LatLngPoly(outer, *holes)


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
        h3_poly = _geojson_to_h3poly(gdf.geometry[0].__geo_interface__)
        hexs = h3.polygon_to_cells(h3_poly, resolution)
        all_polys = gpd.GeoSeries(
            list(map(_hex_to_polygon, hexs)), index=list(hexs), crs="EPSG:4326"
        )

        gdf = gpd.GeoDataFrame(geometry=all_polys, crs="EPSG:4326")
        return gdf

    elif isinstance(gdf.geometry.iloc[0], MultiPolygon):
        parts_lst = []
        for _, row in gdf.explode(index_parts=True).loc[0].iterrows():
            h3_poly = _geojson_to_h3poly(row.geometry.__geo_interface__)
            hexs = h3.polygon_to_cells(h3_poly, resolution)
            all_polys = gpd.GeoSeries(
                list(map(_hex_to_polygon, hexs)), index=list(hexs), crs="EPSG:4326"
            )

            part_gdf = gpd.GeoDataFrame(geometry=all_polys)
            parts_lst.append(part_gdf)

        return pd.concat(parts_lst)
