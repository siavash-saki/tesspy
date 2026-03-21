"""
Unit tests for Voronoi tessellation functions (no network required).
"""

import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon

from tesspy.methods.voronoi import voronoi_polygons


def test_returns_polygons_for_simple_points():
    """voronoi_polygons yields one Polygon per input point."""
    points = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, 0.5],
            [0.0, 1.0],
            [1.0, 1.0],
        ]
    )
    vor = Voronoi(points)
    result = list(voronoi_polygons(vor, diameter=2.0))

    assert len(result) == len(points)
    for poly in result:
        assert isinstance(poly, Polygon)


def test_all_polygons_are_valid():
    """All yielded polygons should be geometrically valid."""
    points = np.array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [1.0, 2.0],
            [1.0, 0.5],
        ]
    )
    vor = Voronoi(points)
    result = list(voronoi_polygons(vor, diameter=5.0))

    for poly in result:
        assert poly.is_valid


def test_polygons_have_nonzero_area():
    """Each Voronoi polygon should have positive area."""
    points = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.5, 0.5],
        ]
    )
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


def test_diameter_covers_wide_generator_span():
    """Voronoi polygons must cover the full bounding box when diameter matches the spread.

    Regression test for tessellation.py hardcoding diameter=0.1.  For generators
    spanning more than 0.1 units, infinite Voronoi regions were truncated too early,
    leaving coverage gaps near the edges of the generator set.

    Verifies that a diameter equal to the bounding-box diagonal produces complete
    coverage, while diameter=0.1 fails for widely-spaced generators.
    """
    from shapely.geometry import box
    from shapely.ops import unary_union

    # Generators spanning 2.0 units — much larger than the old hardcoded 0.1
    points = np.array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [0.0, 2.0],
            [2.0, 2.0],
            [1.0, 1.0],
        ]
    )
    vor = Voronoi(points)
    bounding_box = box(0.0, 0.0, 2.0, 2.0)

    # With the correct diameter (diagonal ≈ 2.83), coverage should be complete.
    x_range = points[:, 0].max() - points[:, 0].min()
    y_range = points[:, 1].max() - points[:, 1].min()
    diameter = np.sqrt(x_range**2 + y_range**2)
    polygons = list(voronoi_polygons(vor, diameter))

    union = unary_union(polygons)
    assert union.covers(bounding_box), (
        "Voronoi union does not cover bounding box — diameter may be too small"
    )

    # With the old hardcoded diameter=0.1, coverage fails for wide generator spread
    small_polygons = list(voronoi_polygons(vor, 0.1))
    small_union = unary_union(small_polygons)
    assert not small_union.covers(bounding_box), (
        "Expected coverage gap with diameter=0.1 for wide generator spread"
    )
