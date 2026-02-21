"""
Data retrieval layer: POI and road network from OpenStreetMap.

This subpackage encapsulates all external network calls so the rest of the
package can stay free of I/O concerns.

Modules
-------
poi      : POIdata class — queries OSM Overpass API for Points of Interest
roads    : RoadData class — downloads road networks via osmnx
_geo     : geocoding helper (get_city_polygon) and count_poi_per_tile utility
_overpass: coordinate rounding helpers for Overpass bbox queries

Exports
-------
POIdata, RoadData, get_city_polygon, count_poi_per_tile

Depends on
----------
requests (Overpass API), osmnx (geocoding + road network), geopandas, shapely,
tesspy._constants (OSM_PRIMARY_FEATURES, OSM_HIGHWAY_TYPES, DEFAULT_POI_CATEGORIES)
"""

from tesspy.data._geo import count_poi_per_tile, get_city_polygon
from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData

__all__ = ["POIdata", "RoadData", "get_city_polygon", "count_poi_per_tile"]
