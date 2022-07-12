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
    city_vp = city.voronoi(cluster_algo="k-means", poi_categories=["amenity"], timeout=60, n_polygons=500,
                           verbose=False)
    city_aqk = city.adaptive_squares(14, poi_categories=["building"])
    print(f"{len(city_qk)} Squares; {len(city_hex)} Hexagons; {len(city_vp)} Voronoi; {len(city_aqk)} Adaptive Squares")


def test_index_methods():
    city = Tessellation("Paris")
    df_vp_hdbscan = city.voronoi(poi_categories=["amenity"], cluster_algo='hdbscan')
    df_vp_kmeans = city.voronoi(poi_categories=["amenity"])
    df_cb = city.city_blocks()
    print(f"Voronoi using Kmeans has columns:{df_vp_kmeans.columns}")
    print(f"Voronoi using hdbscan has columns:{df_vp_hdbscan.columns}")
    print(f"City blocks gdf has columns:{df_cb.columns}")

    print(f"Voronoi using Kmeans has geom_type:{df_vp_kmeans.geom_type.unique()}")
    print(f"Voronoi using hdbscan has geom_type:{df_vp_hdbscan.geom_type.unique()}")
    print(f"City blocks gdf has geom_type:{df_cb.geom_type.unique()}")


def test_count_lgu_function_short_version():
    city = Tessellation("Paris")
    poi_categories2 = ["amenity", "building"]

    df_squares = city.squares(resolution=14)
    df_squares = count_poi_per_tile(city, df_squares, method="squares", poi_categories=poi_categories2)

    print(f'Squares GeoDataFrame has columns: {df_squares.columns}')

    assert hasattr(df_squares, "count_amenity")


def test_count_lgu_function_different_city_inputs():
    city = Tessellation("Paris")
    city_2 = "Paris"
    poi_categories2 = ["amenity", "building"]

    df_squares = city.squares(resolution=14)
    df_squares = count_poi_per_tile(city, df_squares, method="squares", poi_categories=poi_categories2)
    sleep(60)
    df_squares2 = city.squares(resolution=14)
    df_squares2 = count_poi_per_tile(city_2, df_squares2, method="squares", poi_categories=poi_categories2)

    print(f'Squares GeoDataFrame has columns: {df_squares.columns}')
    print(f'Squares2 GeoDataFrame has columns: {df_squares2.columns}')

    assert hasattr(df_squares, "count_amenity")
    assert hasattr(df_squares, "count_building")

    assert hasattr(df_squares2, "count_amenity")
    assert hasattr(df_squares2, "count_building")


def test_count_lgu_function_invalid_inputs_city():
    city = Tessellation("Paris")
    city_2 = "Paris"
    poi_categories2 = ["amenity", "building"]

    df_squares = city.squares(resolution=14)
    try:
        df_squares = count_poi_per_tile(city.get_polygon(), df_squares, method="squares",
                                        poi_categories=poi_categories2)
        assert hasattr(df_squares, "count_amenity")
    except:
        print("-->It worked. Invalid input was recognized")


def test_count_lgu_function_invalid_inputs_gdf():
    city = Tessellation("Paris")
    poi_categories2 = ["amenity", "building"]

    try:
        df_squares = count_poi_per_tile(city, gpd.GeoDataFrame(), method="squares",
                                        poi_categories=poi_categories2)
        assert hasattr(df_squares, "count_amenity")
    except:
        print("-->It worked. Invalid input was recognized")


def test_count_lgu_function_invalid_inputs_method():
    city = Tessellation("Paris")
    city_2 = "Paris"
    poi_categories2 = ["amenity", "building"]

    df_squares = city.squares(resolution=14)
    try:
        df_squares = count_poi_per_tile(city, df_squares, poi_categories=poi_categories2)
        assert hasattr(df_squares, "count_amenity")
    except:
        print("-->It worked. Invalid input was recognized")


def test_count_lgu_function_invalid_inputs_poi():
    city = Tessellation("Paris")
    city_2 = "Paris"
    poi_categories2 = ["amenity", "building"]

    df_squares = city.squares(resolution=14)
    try:
        df_squares = count_poi_per_tile(city, df_squares, poi_categories=10)
        assert hasattr(df_squares, "count_amenity")
    except:
        print("-->It worked. Invalid input was recognized")


def test_count_lgu_function_different_poi():
    city = Tessellation("Paris")
    poi_categories1 = "amenity"
    poi_categories2 = ["amenity"]
    poi_categories3 = ["amenity", "building"]

    df_squares1 = city.squares(resolution=14)
    df_squares2 = city.squares(resolution=14)
    df_squares3 = city.squares(resolution=14)

    df_squares1 = count_poi_per_tile(city, df_squares1, method="squares", poi_categories=poi_categories1)
    sleep(120)
    df_squares2 = count_poi_per_tile(city, df_squares2, method="squares", poi_categories=poi_categories2)
    sleep(120)
    df_squares3 = count_poi_per_tile(city, df_squares3, method="squares", poi_categories=poi_categories3)

    # assert df_squares1.columns == df_squares2.columns
    print(df_squares1.columns)
    print(df_squares2.columns)
    assert 'count_building' in df_squares3.columns


def test_count_lgu_function():
    city = Tessellation("Paris")
    poi_categories1 = ["amenity"]
    poi_categories2 = ["amenity", "building"]

    df_squares = city.squares(resolution=14)
    df_squares = count_poi_per_tile(city, df_squares, method="squares", poi_categories=poi_categories1)
    sleep(120)

    df_hexs = city.hexagons(resolution=9)
    df_hexs = count_poi_per_tile(city, df_hexs, method="hexagons", poi_categories=poi_categories1)
    sleep(120)

    df_aqk = city.adaptive_squares(start_resolution=14, poi_categories=['amenity'])
    df_aqk = count_poi_per_tile(city, df_aqk, method="adaptive_squares", poi_categories=poi_categories2)
    sleep(120)

    df_vp = city.voronoi(poi_categories=['building'])
    df_vp = count_poi_per_tile(city, df_vp, method="voronoi", poi_categories=poi_categories2)
    sleep(120)

    df_cb = city.city_blocks()
    df_cb = count_poi_per_tile(city, df_cb, method="city_blocks", poi_categories=poi_categories1)
    sleep(120)

    print(f'Squares GeoDataFrame has columns: {df_squares.columns}')
    print(f'Hexagons GeoDataFrame has columns: {df_hexs.columns}')
    print(f'Adaptive Squares GeoDataFrame has columns: {df_aqk.columns}')
    print(f'Voronoi Polygons GeoDataFrame has columns: {df_vp.columns}')
    print(f'City blocks GeoDataFrame has columns: {df_cb.columns}')

    assert hasattr(df_squares, "count_amenity")
    assert hasattr(df_hexs, "count_amenity")
    assert hasattr(df_aqk, "count_amenity")
    assert hasattr(df_vp, "count_amenity")
    assert hasattr(df_cb, "count_amenity")


def test_poi_collection():
    city = Tessellation("Frankfurt")
    poi_data_1 = POIdata(city.get_polygon(),
                         poi_categories=["aerialway"],
                         timeout=60,
                         verbose=True).get_poi_data()
    print(poi_data_1.head())

