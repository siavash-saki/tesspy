"""
Clustering utility for city blocks hierarchical merging.
"""

import numpy as np
from sklearn.cluster import AgglomerativeClustering


def _count_clusters(coordinates: np.ndarray, dist_threshold: int) -> int:
    """Fit AgglomerativeClustering and return the number of clusters."""
    model = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=dist_threshold,
        metric="euclidean",
        compute_full_tree=True,
    )
    model.fit(coordinates)
    labels = model.labels_
    return len(set(labels)) - (1 if -1 in labels else 0)


def get_hierarchical_clustering_parameter(
    coordinates: np.ndarray, threshold: int
) -> int | None:
    """
    Find a distance_threshold for AgglomerativeClustering that yields
    fewer clusters than the given threshold.

    Uses binary search over [200, 1200] to minimize the number of model fits.

    Parameters
    ----------
    coordinates : numpy.ndarray
        Array of (x, y) centroid coordinates
    threshold : int
        Maximum acceptable number of clusters

    Returns
    --------
    th : int or None
        The distance_threshold value, or None if no suitable value was found
        in the search range (200–1200)
    """
    low, high = 200, 1200

    # Check that the upper bound is sufficient
    if _count_clusters(coordinates, high) >= threshold:
        return None

    # Check that the lower bound isn't already sufficient
    if _count_clusters(coordinates, low) < threshold:
        return low

    # Binary search for the smallest threshold (in steps of 100)
    while high - low > 100:
        mid = ((low + high) // 200) * 100  # snap to nearest 100
        if mid == low:
            mid = low + 100
        if _count_clusters(coordinates, mid) < threshold:
            high = mid
        else:
            low = mid

    return high
