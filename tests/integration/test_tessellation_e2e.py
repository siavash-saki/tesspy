"""
Integration tests for the Tessellation class (requires OSM API access).
These tests make real network calls and may be slow.
"""

import geopandas as gpd
import pytest

from tesspy.data._geo import count_poi_per_tile, get_city_polygon
from tesspy.tessellation import Tessellation
from tests.conftest import call_with_osm_retry


@pytest.mark.integration
def test_get_city_polygon():
    soho_poly = get_city_polygon("SOHO, London")

    assert isinstance(soho_poly, gpd.GeoDataFrame)
    assert len(soho_poly) == 1
    assert soho_poly.crs is not None

    geom = soho_poly["geometry"].iloc[0]
    assert geom.is_valid
    assert geom.geom_type in ("Polygon", "MultiPolygon")

    # SOHO centroid should be roughly at (-0.135, 51.513)
    centroid = geom.centroid
    assert -0.15 < centroid.x < -0.12
    assert 51.50 < centroid.y < 51.52


@pytest.mark.integration
def test_squares_city_1():
    city1 = Tessellation("Hamburg")
    for i in range(7, 15):
        city1_sq = city1.squares(i)
        city1_sq_2 = city1.squares(i + 1)
        assert len(city1_sq) < len(city1_sq_2)


@pytest.mark.integration
def test_squares_city_2():
    city2 = Tessellation("Key West")
    city2_sq = city2.squares(15)

    assert isinstance(city2_sq, gpd.GeoDataFrame)
    assert len(city2_sq) > 0
    assert hasattr(city2_sq, "geometry")
    assert hasattr(city2_sq, "quadkey")


@pytest.mark.integration
def test_squares_country():
    country = Tessellation("Sri Lanka")
    country_sq = country.squares(12)

    assert isinstance(country_sq, gpd.GeoDataFrame)
    assert len(country_sq) > 0
    assert hasattr(country_sq, "geometry")
    assert hasattr(country_sq, "quadkey")


@pytest.mark.integration
def test_hexagons_city1():
    city1 = Tessellation("Hamburg")
    for i in range(5, 8):
        city1_hex = city1.hexagons(i)
        city2_hex = city1.hexagons(i + 1)
        assert len(city1_hex) < len(city2_hex)


@pytest.mark.integration
def test_hexagons_city2():
    city2 = Tessellation("Key West")
    city2_hx = city2.hexagons(8)

    assert isinstance(city2_hx, gpd.GeoDataFrame)
    assert len(city2_hx) > 0
    assert hasattr(city2_hx, "geometry")
    assert hasattr(city2_hx, "hex_id")


@pytest.mark.integration
def test_hexagons_country():
    country = Tessellation("Sri Lanka")
    country_hx = country.hexagons(4)

    assert isinstance(country_hx, gpd.GeoDataFrame)
    assert len(country_hx) > 0
    assert hasattr(country_hx, "geometry")
    assert hasattr(country_hx, "hex_id")


@pytest.mark.integration
@pytest.mark.slow
def test_adaptive_squares():
    city = Tessellation("Mainz")
    city_sq = city.squares(14)

    city_asq = call_with_osm_retry(
        city.adaptive_squares, 14, poi_categories=["leisure"], timeout=60, verbose=False
    )

    assert isinstance(city_asq, gpd.GeoDataFrame)
    assert len(city_asq) > 0
    assert (city_asq["count"] < 1).sum() == 0
    assert hasattr(city_asq, "geometry")
    assert hasattr(city_asq, "quadkey")
    assert hasattr(city_asq, "count")
    assert len(city_asq) > len(city_sq)


@pytest.mark.integration
@pytest.mark.slow
def test_voronoi():
    city = Tessellation("Mainz")

    city_km = call_with_osm_retry(
        city.voronoi,
        cluster_algo="k-means",
        poi_categories=["leisure"],
        timeout=60,
        n_polygons=20,
        verbose=False,
    )

    city_h = city.voronoi(
        cluster_algo="hdbscan",
        poi_categories=["leisure"],
        timeout=60,
        n_polygons=20,
        verbose=False,
    )

    assert isinstance(city_km, gpd.GeoDataFrame)
    assert len(city_km) > 0
    assert hasattr(city_km, "voronoi_id")

    assert isinstance(city_h, gpd.GeoDataFrame)
    assert len(city_h) > 0
    assert hasattr(city_h, "voronoi_id")


@pytest.mark.integration
@pytest.mark.slow
def test_city_blocks():
    city = Tessellation("Innenstadt, Frankfurt")
    city_cb = city.city_blocks(n_polygons=100, verbose=False)

    assert isinstance(city_cb, gpd.GeoDataFrame)
    assert len(city_cb) > 0
    assert hasattr(city_cb, "geometry")
    assert hasattr(city_cb, "cityblock_id")


@pytest.mark.integration
@pytest.mark.slow
def test_count_poi_per_tile_basic():
    city = Tessellation("Nizza")
    df_squares = city.squares(resolution=14)

    df_with_counts = call_with_osm_retry(
        count_poi_per_tile, "Nizza", df_squares, poi_categories=["leisure"]
    )

    assert isinstance(df_with_counts, gpd.GeoDataFrame)
    assert len(df_with_counts) > 0
    assert hasattr(df_with_counts, "leisure")
    assert df_with_counts["leisure"].min() == 0
    assert df_with_counts["leisure"].max() > 0


def test_count_poi_per_tile_invalid_city_type():
    """Passing a non-string as city should raise ValueError immediately (no network)."""
    from shapely.geometry import Polygon

    poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    fake_gdf = gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")

    with pytest.raises(ValueError):
        count_poi_per_tile(fake_gdf, fake_gdf)


@pytest.mark.integration
def test_count_poi_per_tile_empty_gdf():
    """Passing an empty GeoDataFrame should raise ValueError."""
    with pytest.raises(ValueError):
        count_poi_per_tile("Nizza", gpd.GeoDataFrame())


@pytest.mark.integration
def test_count_poi_per_tile_invalid_poi_type():
    """Passing an int as poi_categories should raise ValueError."""
    with pytest.raises(ValueError):
        count_poi_per_tile("Nizza", gpd.GeoDataFrame(geometry=[]), poi_categories=10)
