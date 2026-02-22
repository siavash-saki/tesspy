"""
Tessellation algorithm implementations.

Each module contains the core functions for one tessellation method.
These are pure computational functions that operate on GeoDataFrames and do
not make network requests.  The Tessellation class in tesspy.tessellation
orchestrates them.

Modules
-------
squares     : regular and adaptive square grids (mercantile quadkeys)
hexagons    : H3 hexagon grids (Uber H3)
voronoi     : Voronoi polygon construction (scipy)
city_blocks : road-network-based city block polygons (osmnx + shapely)
_clustering : hierarchical clustering helper (scikit-learn)

Exports
-------
get_squares_polyfill, get_adaptive_squares, count_poi,
get_h3_hexagons, voronoi_polygons,
create_city_blocks, merge_city_blocks,
split_linestring, explode, create_blocks, get_rest_polygon (deprecated),
get_hierarchical_clustering_parameter

Depends on
----------
geopandas, shapely, mercantile, h3, scipy, scikit-learn
"""

from tesspy.methods._clustering import get_hierarchical_clustering_parameter
from tesspy.methods.city_blocks import (
    create_blocks,
    create_city_blocks,
    explode,
    get_rest_polygon,
    merge_city_blocks,
    split_linestring,
)
from tesspy.methods.hexagons import get_h3_hexagons
from tesspy.methods.squares import count_poi, get_adaptive_squares, get_squares_polyfill
from tesspy.methods.voronoi import voronoi_polygons

__all__ = [
    "count_poi",
    "get_squares_polyfill",
    "get_adaptive_squares",
    "get_h3_hexagons",
    "voronoi_polygons",
    "create_city_blocks",
    "merge_city_blocks",
    "split_linestring",
    "explode",
    "create_blocks",
    "get_rest_polygon",
    "get_hierarchical_clustering_parameter",
]
