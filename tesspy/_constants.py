"""
Shared OSM constants used across the package.
These lists serve as the single source of truth to avoid duplication.

Exports
-------
OSM_PRIMARY_FEATURES : list of all 29 OSM primary map feature tag keys
OSM_HIGHWAY_TYPES : list of 19 OSM highway type values, ordered by road class
DEFAULT_POI_CATEGORIES : default POI categories used when none are specified
"""

# Default POI categories used when no explicit poi_categories argument is given.
# Defined here once so every method that has a poi_categories parameter
# refers to the same value rather than a copy-pasted literal.
DEFAULT_POI_CATEGORIES: list[str] = ["amenity", "building"]

OSM_PRIMARY_FEATURES = [
    "aerialway",
    "aeroway",
    "amenity",
    "barrier",
    "boundary",
    "building",
    "craft",
    "emergency",
    "geological",
    "healthcare",
    "highway",
    "historic",
    "landuse",
    "leisure",
    "man_made",
    "military",
    "natural",
    "office",
    "place",
    "power",
    "public_transport",
    "railway",
    "route",
    "shop",
    "sport",
    "telecom",
    "tourism",
    "water",
    "waterway",
]

OSM_HIGHWAY_TYPES = [
    "motorway",
    "trunk",
    "primary",
    "secondary",
    "tertiary",
    "residential",
    "unclassified",
    "motorway_link",
    "trunk_link",
    "primary_link",
    "secondary_link",
    "living_street",
    "pedestrian",
    "track",
    "bus_guideway",
    "footway",
    "path",
    "service",
    "cycleway",
]
