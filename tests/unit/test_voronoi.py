"""
Unit tests for Voronoi tessellation functions (no network required).
"""

import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon

from tesspy.methods.voronoi import voronoi_polygons


def test_returns_polygons_for_simple_points():
    """voronoi_polygons yields one Polygon per input point."""
    points = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.5, 1.0],
        [0.5, 0.5],
        [0.0, 1.0],
        [1.0, 1.0],
    ])
    vor = Voronoi(points)
    result = list(voronoi_polygons(vor, diameter=2.0))

    assert len(result) == len(points)
    for poly in result:
        assert isinstance(poly, Polygon)


def test_all_polygons_are_valid():
    """All yielded polygons should be geometrically valid."""
    points = np.array([
        [0.0, 0.0],
        [2.0, 0.0],
        [1.0, 2.0],
        [1.0, 0.5],
    ])
    vor = Voronoi(points)
    result = list(voronoi_polygons(vor, diameter=5.0))

    for poly in result:
        assert poly.is_valid


def test_polygons_have_nonzero_area():
    """Each Voronoi polygon should have positive area."""
    points = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [0.5, 0.5],
    ])
    vor = Voronoi(points)
    result = list(voronoi_polygons(vor, diameter=5.0))

    for poly in result:
        assert poly.area > 0


def test_is_generator():
    """voronoi_polygons should be a generator, not return a list."""
    points = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])
    vor = Voronoi(points)
    gen = voronoi_polygons(vor, diameter=2.0)

    # Should be a generator, consumable with next()
    first = next(gen)
    assert isinstance(first, Polygon)
