from tesspy.tessellation import *
from tesspy.tessellation_functions import *
from tesspy.poi_data import *
import geopandas as gpd
from shapely import wkt
from time import sleep


def test_get_city_polygon():
    soho_boundary = wkt.loads(
        "POLYGON ((-0.1418295 51.5150964, -0.1416219 51.5144947, -0.1412154 51.513919, -0.1407406 51.5133165, -0.1401938 51.5127226, -0.1396956 51.5121827, "
        "-0.139218 51.5116136, -0.1387082 51.511047, -0.1381347 51.510443, -0.138056 51.5103684, -0.1379041 51.5102851, -0.1376615 51.5101849, "
        "-0.1374938 51.5101273, -0.1373739 51.5100897, -0.137174 51.5100423, -0.1369607 51.5100018, -0.1367699 51.5099796, -0.1365341 51.5099691, "
        "-0.1363282 51.5099663, -0.136199 51.5099763, -0.1359833 51.5099983, -0.1358617 51.5100121, -0.1357304 51.5100291, -0.1354012 51.5101019, "
        "-0.1352774 51.510131, -0.1352627 51.5101223, -0.13512 51.5101604, -0.1350045 51.5101994, -0.134819 51.5102848, -0.1345553 51.5102187, "
        "-0.1344994 51.510216, -0.1344519 51.51022, -0.1344024 51.5102344, -0.1343994 51.5102367, -0.1343726 51.5102578, -0.1343369 51.5103096, "
        "-0.1342521 51.510441, -0.1340647 51.5106759, -0.1339251 51.5108473, -0.1338178 51.5109648, -0.133485 51.511293, -0.1331837 51.5114941, "
        "-0.1330286 51.5116095, -0.1329784 51.5116414, -0.1328557 51.5116994, -0.1325545 51.5118259, -0.1322504 51.5119561, -0.1314392 51.5122874, "
        "-0.1304367 51.5126664, -0.1299327 51.5128501, -0.1297379 51.5129546, -0.1294894 51.5130997, -0.1293784 51.5131726, -0.1294407 51.5131871, "
        "-0.1294805 51.5132494, -0.1295057 51.5132957, -0.1295108 51.5133475, -0.1294837 51.5134495, -0.1294596 51.5134934, -0.1294072 51.5135236, "
        "-0.1294217 51.5135966, -0.1294362 51.5136317, -0.1295306 51.5137946, -0.129567 51.5138695, -0.129636 51.5139849, -0.1296451 51.5140012, "
        "-0.1296867 51.5141046, -0.1297291 51.5142121, -0.1297884 51.5143634, -0.1298403 51.5145007, -0.1299175 51.514682, -0.1300152 51.51491, "
        "-0.1301523 51.5152675, -0.1306884 51.5151383, -0.1307937 51.5153197, -0.1308895 51.5154833, -0.1307956 51.5155047, -0.1308833 51.5155937, "
        "-0.1309068 51.5155977, -0.1309316 51.5155987, -0.1309506 51.5155934, -0.1311246 51.5158882, -0.1312572 51.5160966, -0.1312752 51.5161264, "
        "-0.1309837 51.5161582, -0.1307013 51.5161832, -0.1307307 51.516365, -0.1313097 51.5163318, -0.1321982 51.5162733, -0.1329511 51.516213, "
        "-0.1335465 51.5161508, -0.1353119 51.5159597, -0.1360908 51.5158808, -0.1376583 51.5157261, -0.1386282 51.5156109, -0.1393623 51.5155356, "
        "-0.1401402 51.5154556, -0.1411432 51.5153355, -0.1415848 51.5152797, -0.1417661 51.5152281, -0.1418295 51.5150964))"
    )

    soho_poly = get_city_polygon("SOHO, London")
    soho_poly = soho_poly["geometry"].iloc[0]

    assert soho_poly == soho_boundary


def test_squares_city_1():
    # different resolutions -> increasing resolution produces more squares
    city1 = Tessellation("Hamburg")
    for i in range(7, 15):
        city1_sq = city1.squares(i)
        city1_sq_2 = city1.squares(i + 1)
        assert len(city1_sq) < len(city1_sq_2)


def test_squares_city_2():
    # MultiPolygon like Islands
    city2 = Tessellation("Key West")
    city2_sq = city2.squares(15)

    assert type(city2_sq) == gpd.GeoDataFrame
    assert len(city2_sq) > 0
    assert hasattr(city2_sq, "geometry")
    assert hasattr(city2_sq, "quadkey")


def test_squares_country():
    # Whole country
    country = Tessellation("Sri Lanka")
    country_sq = country.squares(12)

    assert type(country_sq) == gpd.GeoDataFrame
    assert len(country_sq) > 0
    assert hasattr(country_sq, "geometry")
    assert hasattr(country_sq, "quadkey")


def test_hexagons_city1():
    # different resolutions -> increasing resolution produces more squares
    city1 = Tessellation("Hamburg")
    for i in range(5, 8):
        city1_hex = city1.hexagons(i)
        city2_hex = city1.hexagons(i + 1)
        assert len(city1_hex) < len(city2_hex)


def test_hexagons_city2():
    # MultiPolygon like Islands
    city2 = Tessellation("Key West")
    city2_hx = city2.hexagons(8)

    assert type(city2_hx) == gpd.GeoDataFrame
    assert len(city2_hx) > 0
    assert hasattr(city2_hx, "geometry")
    assert hasattr(city2_hx, "hex_id")


def test_hexagons_country():
    # Whole country
    country = Tessellation("Sri Lanka")
    country_hx = country.hexagons(4)

    assert type(country_hx) == gpd.GeoDataFrame
    assert len(country_hx) > 0
    assert hasattr(country_hx, "geometry")
    assert hasattr(country_hx, "hex_id")


def test_adaptive_squares():
    city = Tessellation("Mainz")
    city_sq = city.squares(14)

    while True:
        try:
            city_asq = city.adaptive_squares(
                14, poi_categories=["leisure"], timeout=60, verbose=False
            )
            break
        except RuntimeError:
            sleep(5 * 60)
            continue

    assert type(city_asq) == gpd.GeoDataFrame
    assert len(city_asq) > 0
    assert (city_asq["count"] < 1).sum() == 0
    assert hasattr(city_asq, "geometry")
    assert hasattr(city_asq, "quadkey")
    assert hasattr(city_asq, "count")
    assert len(city_asq) > len(city_sq)


def test_voronoi():
    city = Tessellation("Mainz")

    while True:
        try:
            city_km = city.voronoi(
                cluster_algo="k-means",
                poi_categories=["leisure"],
                timeout=60,
                n_polygons=20,
                verbose=False,
            )
            break
        except RuntimeError:
            sleep(5 * 60)
            continue

    city_h = city.voronoi(
        cluster_algo="hdbscan",
        poi_categories=["leisure"],
        timeout=60,
        n_polygons=20,
        verbose=False,
    )

    assert type(city_km) == gpd.GeoDataFrame
    assert len(city_km) > 0
    assert hasattr(city_km, "geometry")
    assert hasattr(city_km, "voronoi_id")

    assert type(city_h) == gpd.GeoDataFrame
    assert len(city_h) > 0
    assert hasattr(city_h, "geometry")
    assert hasattr(city_h, "voronoi_id")


def test_city_blocks():
    city = Tessellation("Innenstadt, Frankfurt")
    nb_poly = 100
    sleep(30)
    city_cb = city.city_blocks(n_polygons=nb_poly, verbose=False)

    assert type(city_cb) == gpd.GeoDataFrame
    assert len(city_cb) > 0
    assert hasattr(city_cb, "geometry")
    assert hasattr(city_cb, "cityblock_id")


def test_count_lgu_function_short_version():
    city = Tessellation("Nizza")
    poi_categories2 = ["leisure"]

    df_squares = city.squares(resolution=14)
    while True:
        try:
            df_squares_with_counts = count_poi_per_tile(
                "Nizza", df_squares, poi_categories=poi_categories2
            )
            break
        except RuntimeError:
            sleep(5 * 60)
            continue

    assert type(df_squares_with_counts) == gpd.GeoDataFrame
    assert len(df_squares_with_counts) > 0
    assert hasattr(df_squares_with_counts, "geometry")
    assert hasattr(df_squares_with_counts, "quadkey")
    assert hasattr(df_squares_with_counts, "leisure")
    assert df_squares_with_counts["leisure"].min() == 0
    assert df_squares_with_counts["leisure"].max() > 0
    assert df_squares_with_counts["leisure"].sum() > 0


def test_count_lgu_function_invalid_inputs_city():
    city = Tessellation("Nizza")
    poi_categories2 = ["leisure"]

    df_squares = city.squares(resolution=14)
    try:
        df_squares = count_poi_per_tile(
            city.get_polygon(), df_squares, poi_categories=poi_categories2
        )
        assert hasattr(df_squares, "amenity")
    except ValueError as e:
        assert isinstance(e, ValueError)


def test_count_lgu_function_invalid_inputs_gdf():
    city = Tessellation("Nizza")
    poi_categories2 = ["leisure"]

    try:
        df_squares = count_poi_per_tile(
            city, gpd.GeoDataFrame(), poi_categories=poi_categories2
        )
        assert hasattr(df_squares, "leisure")
    except ValueError as e:
        assert isinstance(e, ValueError)


def test_count_lgu_function_invalid_inputs_method():
    city = Tessellation("Nizza")
    poi_categories2 = ["leisure"]
    df_squares = city.squares(resolution=14)

    try:
        df_squares = count_poi_per_tile(
            city, df_squares, poi_categories=poi_categories2
        )
        assert hasattr(df_squares, "leisure")
    except ValueError as e:
        assert isinstance(e, ValueError)


def test_count_lgu_function_invalid_inputs_poi():
    city = Tessellation("Nizza")
    df_squares = city.squares(resolution=14)

    try:
        df_squares = count_poi_per_tile(city, df_squares, poi_categories=10)
        assert hasattr(df_squares, "amenity")
    except ValueError as e:
        assert isinstance(e, ValueError)


def test_count_lgu_function_different_poi():
    city = Tessellation("Mitte, Berlin")
    poi_categories1 = "telecom"
    poi_categories3 = ["telecom", "public_transport"]
    df_squares = city.squares(resolution=14)
    sleep(180)

    while True:
        try:
            df_squares1 = count_poi_per_tile(
                "Mitte, Berlin", df_squares, poi_categories=poi_categories1
            )
            break
        except RuntimeError:
            sleep(5 * 60)
            continue

    while True:
        try:
            df_squares2 = count_poi_per_tile(
                "Mitte, Berlin", df_squares, poi_categories=poi_categories3
            )
            break
        except RuntimeError:
            sleep(5 * 60)
            continue

    assert hasattr(df_squares1, "telecom")
    assert hasattr(df_squares2, "telecom")
    assert hasattr(df_squares2, "public_transport")

    assert df_squares1["telecom"].sum() > 0
    assert df_squares2["telecom"].sum() > 0
