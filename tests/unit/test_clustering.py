"""
Unit tests for hierarchical clustering parameter search (no network required).
"""

import numpy as np

from tesspy.methods._clustering import get_hierarchical_clustering_parameter


def test_returns_int_for_easy_threshold():
    """With a high threshold, a low distance_threshold should suffice."""
    # 3 tight clusters, well separated
    coords = np.array([
        [0, 0], [1, 1], [2, 2],
        [500, 500], [501, 501], [502, 502],
        [1000, 1000], [1001, 1001], [1002, 1002],
    ])
    result = get_hierarchical_clustering_parameter(coords, threshold=10)

    assert isinstance(result, int)
    assert result >= 200


def test_returns_none_when_impossible():
    """When points are too spread out and threshold is 1, no value can satisfy."""
    # Many well-separated points — even distance_threshold=1200 yields >1 cluster
    rng = np.random.default_rng(42)
    coords = rng.uniform(0, 100_000, size=(200, 2))
    result = get_hierarchical_clustering_parameter(coords, threshold=1)

    assert result is None


def test_returns_first_satisfying_threshold():
    """The function should return the smallest distance_threshold that works."""
    # Two tight clusters, separated by ~600 units
    coords = np.array([
        [0, 0], [1, 0], [0, 1],
        [600, 600], [601, 600], [600, 601],
    ])
    result = get_hierarchical_clustering_parameter(coords, threshold=3)

    assert isinstance(result, int)
    # The first threshold that merges 2 clusters into fewer than 3
    # should be the smallest in the search range that bridges the gap
    assert result in [i * 100 for i in range(2, 13)]


def test_search_range():
    """The function searches thresholds [200, 300, ..., 1200]."""
    # Single tight cluster — even distance_threshold=200 should give 1 cluster < 5
    coords = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])
    result = get_hierarchical_clustering_parameter(coords, threshold=5)

    assert result == 200  # first in search range
