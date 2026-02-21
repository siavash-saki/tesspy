"""
Backward-compatibility shim for tesspy.poi_data.

Classes have been moved to tesspy.data.*
Import from there directly in new code:

    from tesspy.data.poi import POIdata
    from tesspy.data.roads import RoadData
"""

from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData
from tesspy.data._overpass import geom_ceil, geom_floor

__all__ = ["POIdata", "RoadData", "geom_ceil", "geom_floor"]
