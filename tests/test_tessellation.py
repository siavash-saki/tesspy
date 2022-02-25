from tesspy.tessellation import *
from tesspy.tessellation_functions import *
from tesspy.poi_data import *

from time import sleep


def test_squares():
    # different resolutions -> increasing resolution produces more squares
    city1 = Tessellation("Hamburg")
    for i in range(7, 15):
        city1_sq = city1.squares(i)
        city2_sq = city1.squares(i + 1)
        assert len(city1_sq) < len(city2_sq)

    # MultiPolygon like Islands
    city2 = Tessellation("Key West")
    assert hasattr(city2.squares(15), "geometry")

    # Whole country
    city3 = Tessellation("Sri Lanka")
    assert hasattr(city3.squares(6), "geometry")


def test_hexagons():
    # different resolutions -> increasing resolution produces more squares
    city1 = Tessellation("Hamburg")
    for i in range(5, 10):
        city1_hex = city1.hexagons(i)
        city2_hex = city1.hexagons(i + 1)
        assert len(city1_hex) < len(city2_hex)

    # MultiPolygon like Islands
    city2 = Tessellation("Key West")
    assert hasattr(city2.hexagons(8), "geometry")

    # Whole country
    city3 = Tessellation("Sri Lanka")
    assert hasattr(city3.hexagons(4), "geometry")


def test_adaptive_squares():
    city = Tessellation("Mainz")
    city_sq = city.squares(14)
    city_asq1 = city.adaptive_squares(14, poi_categories=["leisure"], timeout=60, verbose=False)
    city_asq2 = city.adaptive_squares(14, poi_categories=["leisure"], threshold=10000000, timeout=60, verbose=False)

    assert len(city_asq1) > len(city_sq)
    assert hasattr(city_asq1, "geometry")
    assert len(city_asq2) == len(city_sq)

    city2 = Tessellation("Sydney")
    city2_asq = city2.adaptive_squares(14, poi_categories=["amenity"], timeout=60, verbose=False)
    assert hasattr(city2_asq, "geometry")

    city3 = Tessellation("Hanau")
    city3_asq = city3.adaptive_squares(14, poi_categories=["amenity", "building", "leisure"], timeout=60, verbose=False)
    assert hasattr(city3_asq, "geometry")


def test_voronoi():
    city = Tessellation("Buenos Aires")
    city_km = city.voronoi(cluster_algo="k-means",
                           poi_categories=["amenity"],
                           timeout=60,
                           n_polygons=300,
                           verbose=False)
    assert len(city_km) <= 300

    sleep(60)

    city = Tessellation("Buenos Aires")
    city_hdb = city.voronoi(cluster_algo="hdbscan",
                            poi_categories=["building"],
                            timeout=60,
                            min_cluster_size=10,
                            verbose=False)
    assert hasattr(city_hdb, "geometry")


def test_city_blocks():
    city = Tessellation("Frankfurt")
    nb_LGUs = 1800
    city_cb = city.city_blocks(n_polygons=nb_LGUs, verbose=True)
    print(f"We have {len(city_cb)} - the userinput was {nb_LGUs}")
    assert hasattr(city_cb, "geometry")


def test_new_shapely():
    city = Tessellation("Key West")
    city_qk = city.squares(15)
    city_hex = city.hexagons(8)
    city_vp = city.voronoi(cluster_algo="k-means",poi_categories=["amenity"], timeout=60, n_polygons=500, verbose=False)
    city_aqk = city.adaptive_squares(14, poi_categories=["building"])
    print(f"{len(city_qk)} Squares; {len(city_hex)} Hexagons; {len(city_vp)} Voronoi; {len(city_aqk)} Adaptive Squares")
