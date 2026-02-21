"""
Internal helpers for Overpass API coordinate rounding.
"""

import numpy as np


def geom_ceil(coordinate: float, precision: int = 4) -> float:
    """Round coordinate up to specified decimal precision."""
    return float(np.true_divide(np.ceil(coordinate * 10**precision), 10**precision))


def geom_floor(coordinate: float, precision: int = 4) -> float:
    """Round coordinate down to specified decimal precision."""
    return float(np.true_divide(np.floor(coordinate * 10**precision), 10**precision))
