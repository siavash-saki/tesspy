"""
Backward-compatibility shim for tesspy.tessellation_functions.

.. deprecated::
    This module is deprecated and will be removed in v0.3.0.
    Import from tesspy.methods directly:

    from tesspy.methods.squares import (
        count_poi, get_squares_polyfill, get_adaptive_squares,
    )
    from tesspy.methods.hexagons import get_h3_hexagons
    from tesspy.methods.voronoi import voronoi_polygons
    from tesspy.methods.city_blocks import (
        split_linestring, explode, create_blocks, get_rest_polygon,
    )
    from tesspy.methods._clustering import get_hierarchical_clustering_parameter
"""

import warnings

warnings.warn(
    "tesspy.tessellation_functions is deprecated and will be removed in v0.3.0. "
    "Import from tesspy.methods instead: "
    "'from tesspy.methods.squares import get_squares_polyfill, ...'",
    DeprecationWarning,
    stacklevel=2,
)

from tesspy.methods._clustering import (  # noqa: E402
    get_hierarchical_clustering_parameter,
)
from tesspy.methods.city_blocks import (  # noqa: E402
    create_blocks,
    explode,
    get_rest_polygon,
    split_linestring,
)
from tesspy.methods.hexagons import get_h3_hexagons  # noqa: E402
from tesspy.methods.squares import (  # noqa: E402
    count_poi,
    get_adaptive_squares,
    get_squares_polyfill,
)
from tesspy.methods.voronoi import voronoi_polygons  # noqa: E402

__all__ = [
    "count_poi",
    "get_squares_polyfill",
    "get_adaptive_squares",
    "get_h3_hexagons",
    "voronoi_polygons",
    "split_linestring",
    "explode",
    "create_blocks",
    "get_rest_polygon",
    "get_hierarchical_clustering_parameter",
]
