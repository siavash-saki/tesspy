"""
Backward-compatibility shim for tesspy.poi_data.

.. deprecated::
    This module is deprecated and will be removed in v0.3.0.
    Import from tesspy.data directly:

    from tesspy.data.poi import POIdata
    from tesspy.data.roads import RoadData
"""

import warnings

warnings.warn(
    "tesspy.poi_data is deprecated and will be removed in v0.3.0. "
    "Import from tesspy.data instead: "
    "'from tesspy.data.poi import POIdata'",
    DeprecationWarning,
    stacklevel=2,
)

from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData
from tesspy.data._overpass import geom_ceil, geom_floor

__all__ = ["POIdata", "RoadData", "geom_ceil", "geom_floor"]
