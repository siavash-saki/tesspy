"""
Input validation helpers for the Tessellation class.
"""

import logging

import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon

from tesspy.methods.city_blocks import explode

logger = logging.getLogger(__name__)


def _check_input_geodataframe(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Validate that a GeoDataFrame is suitable as a study area input.

    Checks:
    - Exactly one geometry element
    - geometry column is present
    - Geometry is Polygon or MultiPolygon
    - CRS is defined; converts to EPSG:4326 if needed

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        Validated (and possibly reprojected) GeoDataFrame
    """
    if len(gdf) != 1:
        raise ValueError("GeoDataFrame must have exactly one geometry element")

    if not hasattr(gdf, "geometry"):
        raise TypeError("Geometry column missing in GeoDataFrame")

    if type(gdf["geometry"].iloc[0]) not in [Polygon, MultiPolygon]:
        raise TypeError("Geometry must be of type shapely Polygon or MultiPolygon")

    if gdf.crs is None:
        raise ValueError("GeoDataFrame must have a CRS")

    if gdf.crs != "epsg:4326":
        return gdf.to_crs(epsg=4326)

    return gdf


def _check_valid_geometry_gdf(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Ensure all geometries in the GeoDataFrame are single Polygons.
    Any MultiPolygons are exploded into individual Polygon rows.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame

    Returns
    --------
    gdf : geopandas.GeoDataFrame
    """
    if len(gdf) < 1:
        raise ValueError("GeoDataFrame must have at least one geometry element")

    if not hasattr(gdf, "geometry"):
        raise TypeError("Geometry column missing in GeoDataFrame")

    if "MultiPolygon" in gdf.geom_type.unique():
        logger.info("event=geometry.multipolygon.explode")
        gdf = explode(gdf)
        gdf = gdf.reset_index()
        gdf.drop(columns=["level_0", "level_1"], inplace=True)

    return gdf
