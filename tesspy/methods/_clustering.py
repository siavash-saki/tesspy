"""
Clustering utility for city blocks hierarchical merging.
"""

import numpy as np
from sklearn.cluster import AgglomerativeClustering


def get_hierarchical_clustering_parameter(
    coordinates: np.ndarray, threshold: int
) -> int | None:
    """
    Find a distance_threshold for AgglomerativeClustering that yields
    fewer clusters than the given threshold.

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
    dist_threshold = [i * 100 for i in range(2, 13)]
    for th in dist_threshold:
        model = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=th,
            affinity="euclidean",
            compute_full_tree=True,
        )
        model.fit(coordinates)
        labels = model.labels_
        nb_clusters = len(set(labels)) - (1 if -1 in labels else 0)

        if nb_clusters < threshold:
            return th

    return None
