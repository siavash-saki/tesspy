"""
Geocoding helpers and the count_poi_per_tile public utility.
"""

import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point


def get_city_polygon(city: str) -> gpd.GeoDataFrame:
    """
    Retrieve the boundary polygon of a city or region from OSM.

    Parameters
    ----------
    city : str
        Name of a city or address of a region

    Returns
    --------
    df_city : geopandas.GeoDataFrame
        GeoDataFrame containing the boundary polygon
    """
    import osmnx as ox

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df_city = ox.geocode_to_gdf(city)
    df_city = df_city[["osm_id", "geometry"]]
    df_city = df_city.rename(columns={"osm_id": "osmid"})
    return df_city


def count_poi_per_tile(
    city,
    gdf: gpd.GeoDataFrame,
    poi_categories: list[str] | str | None = None,
    timeout: int = 120,
) -> gpd.GeoDataFrame:
    """
    Count POI categories per tessellation tile.

    For each POI category an additional count column is added to the
    tessellation GeoDataFrame. Accepts either a Tessellation object or
    a city name string.

    Parameters
    ----------
    city : tesspy.Tessellation or str
        Tessellation object or city name string for the study area
    gdf : geopandas.GeoDataFrame
        Tessellation GeoDataFrame (output of any Tessellation method)
    poi_categories : list of str or str, default=["amenity", "building"]
        OSM primary map feature categories to count per tile.
        Pass 'all' for all available categories.
    timeout : int, default=120
        TCP timeout in seconds for the OSM Overpass request

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        Tessellation GeoDataFrame with additional count columns per POI category
    """
    # Avoid circular import: Tessellation imports from data, so import here
    from tesspy.tessellation import Tessellation

    if poi_categories is None:
        poi_categories = ["amenity", "building"]

    if type(city) == str:
        city = Tessellation(city)
    else:
        raise ValueError(
            "Please insert a valid city type. Valid types are: "
            "tesspy.Tessellation object or string"
        )

    if len(gdf) < 1:
        raise ValueError(
            "Please insert a valid tessellation GeoDataFrame with at least one tile."
        )

    if type(poi_categories) == str:
        poi_categories = [poi_categories]
    elif type(poi_categories) not in (list,):
        raise ValueError(
            "poi_categories must be a string or list of OSM primary feature names."
        )

    from tesspy.data.poi import POIdata

    df_poi = POIdata(
        city.get_polygon(),
        poi_categories=poi_categories,
        timeout=timeout,
        verbose=False,
    ).get_poi_data()

    points_geom = df_poi[["center_longitude", "center_latitude"]].apply(
        lambda p: Point(p["center_longitude"], p["center_latitude"]), axis=1
    )

    tess_data = gpd.GeoDataFrame(
        geometry=points_geom, data=df_poi[poi_categories], crs="EPSG:4326"
    )

    tess_data["value"] = (
        tess_data.drop(columns=["geometry"]).idxmax(1).where(tess_data.any(1))
    )
    tess_data = tess_data[["value", "geometry"]]

    try:
        idx = [s for s in gdf.columns if s.__contains__("id")][0]
    except IndexError:
        idx = [s for s in gdf.columns if s.__contains__("key")][0]

    spatial_join = gpd.sjoin(gdf, tess_data)
    pivot_table = pd.pivot_table(
        spatial_join, index=idx, columns="value", aggfunc={"value": len}
    )

    pivot_table.columns = pivot_table.columns.droplevel()

    merged_polygons = gdf.merge(pivot_table, how="left", on=idx)
    merged_polygons.fillna(0, inplace=True)

    return merged_polygons
