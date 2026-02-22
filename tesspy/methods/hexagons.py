"""
Hexagon tessellation functions (Uber H3).
"""

import geopandas as gpd
import h3
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
    geom = gdf.geometry.iloc[0]

    if isinstance(geom, Polygon):
        parts = [geom]
    elif isinstance(geom, MultiPolygon):
        parts = list(geom.geoms)
    else:
        raise TypeError(f"Unsupported geometry type: {type(geom)}")

    all_hexs: list[str] = []
    for part in parts:
        h3_poly = _geojson_to_h3poly(part.__geo_interface__)
        all_hexs.extend(h3.polygon_to_cells(h3_poly, resolution))

    # Deduplicate (overlapping parts may share boundary hexagons)
    unique_hexs = list(dict.fromkeys(all_hexs))
    polys = [_hex_to_polygon(h) for h in unique_hexs]

    return gpd.GeoDataFrame(
        geometry=gpd.GeoSeries(polys, index=unique_hexs, crs="EPSG:4326")
    )
