"""
tesspy — geographic tessellation for urban spatial analysis.

tesspy discretises a geographic area into non-overlapping, gap-free subregions
using five tessellation methods driven by OpenStreetMap data:

    * squares         — regular quadkey-based square grid
    * hexagons        — regular H3 hexagon grid (Uber)
    * adaptive_squares — POI-density-driven recursive square subdivision
    * voronoi         — Voronoi polygons seeded by POI cluster centroids
    * city_blocks     — urban blocks formed by the OSM road network

Public API
----------
Tessellation : main class; instantiate with a city name or GeoDataFrame
count_poi_per_tile : count OSM POI categories per tessellation tile
configure_logging : configure package logger output/level

Subpackages
-----------
tesspy.methods : individual tessellation algorithm implementations
tesspy.data    : OSM data retrieval (POI via Overpass API, roads via osmnx)

Example
-------
>>> from tesspy import Tessellation
>>> t = Tessellation("Frankfurt am Main")
>>> squares = t.squares(resolution=12)
>>> hexagons = t.hexagons(resolution=8)
"""

import logging

from tesspy._logging import configure_logging
from tesspy._version import __version__
from tesspy.data._geo import count_poi_per_tile
from tesspy.tessellation import Tessellation

logging.getLogger("tesspy").addHandler(logging.NullHandler())

name = "tesspy"
__author__ = "Siavash Saki and Jonas Hamann"
__author_email__ = "jonas.hamann@fb3.fra-uas.de"

__all__ = [
    "Tessellation",
    "count_poi_per_tile",
    "configure_logging",
    "__version__",
]
