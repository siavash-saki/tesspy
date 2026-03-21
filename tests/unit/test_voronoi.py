"""
Unit tests for Voronoi tessellation functions (no network required).
"""

import numpy as np
import pytest
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


def test_hdbscan_all_noise_raises_value_error(mocker):
    """voronoi(cluster_algo='hdbscan') raises ValueError when all points are noise.

    HDBSCAN assigns label -1 to noise points. When min_cluster_size is too large
    for the dataset, every point gets label -1 and labels_.max() == -1, which
    previously set n_clusters = 0 and crashed inside scipy.spatial.Voronoi.
    The fix should raise a clear ValueError instead.
    """
    import geopandas as gpd
    from shapely.geometry import Polygon

    from tesspy.tessellation import Tessellation

    # Build a minimal Tessellation backed by a real GeoDataFrame (no network)
    area_poly = Polygon([(-0.1, -0.1), (0.1, -0.1), (0.1, 0.1), (-0.1, 0.1)])
    area_gdf = gpd.GeoDataFrame(geometry=[area_poly], crs="EPSG:4326")
    t = Tessellation(area_gdf)

    # Inject a small POI dataframe so the method reaches the HDBSCAN branch.
    # Include the default POI category columns so _get_missing_poi_categories
    # returns [] and no network call is made.
    import pandas as pd

    t.poi_dataframe = pd.DataFrame(
        {
            "center_longitude": [0.0, 0.01, 0.02],
            "center_latitude": [0.0, 0.01, 0.02],
            "amenity": [True, True, True],
            "building": [False, False, False],
        }
    )

    # Patch HDBSCAN so every point is classified as noise (label = -1)
    mock_clustering = mocker.MagicMock()
    mock_clustering.labels_ = np.array([-1, -1, -1])
    mock_hdbscan_cls = mocker.patch(
        "tesspy.tessellation.HDBSCAN", return_value=mock_clustering
    )
    mock_hdbscan_cls.return_value.fit.return_value = mock_clustering

    with pytest.raises(ValueError, match="HDBSCAN found no clusters"):
        t.voronoi(cluster_algo="hdbscan", min_cluster_size=1000)
