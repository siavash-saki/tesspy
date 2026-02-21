"""
Shared test fixtures and utilities for the tesspy test suite.

Unit tests (tests/unit/) run without any network access.
Integration tests (tests/integration/) make real OSM API calls and are
marked with @pytest.mark.integration.
"""

import time
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest
from shapely import wkt

TESTS_DIR = Path(__file__).parent


# ---------------------------------------------------------------------------
# File-based fixtures (no network required)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def leisure_poi_berlin() -> gpd.GeoDataFrame:
    """Pre-downloaded leisure POI GeoDataFrame for Mitte, Berlin."""
    return gpd.read_file(TESTS_DIR / "leisure_poi_Mitte_Berlin.geojson")


@pytest.fixture(scope="session")
def road_data_lille():
    """Pre-downloaded road network GeoDataFrame for Lille."""
    return pd.read_pickle(TESTS_DIR / "road_data_Lille.pkl")


# ---------------------------------------------------------------------------
# Hardcoded polygon fixtures (no network required)
# ---------------------------------------------------------------------------

_SOHO_LONDON_WKT = (
    "POLYGON ((-0.1418295 51.5150964, -0.1416219 51.5144947, -0.1412154 51.513919, "
    "-0.1407406 51.5133165, -0.1401938 51.5127226, -0.1396956 51.5121827, "
    "-0.139218 51.5116136, -0.1387082 51.511047, -0.1381347 51.510443, "
    "-0.138056 51.5103684, -0.1379041 51.5102851, -0.1376615 51.5101849, "
    "-0.1374938 51.5101273, -0.1373739 51.5100897, -0.137174 51.5100423, "
    "-0.1369607 51.5100018, -0.1367699 51.5099796, -0.1365341 51.5099691, "
    "-0.1363282 51.5099663, -0.136199 51.5099763, -0.1359833 51.5099983, "
    "-0.1358617 51.5100121, -0.1357304 51.5100291, -0.1354012 51.5101019, "
    "-0.1352774 51.510131, -0.1352627 51.5101223, -0.13512 51.5101604, "
    "-0.1350045 51.5101994, -0.134819 51.5102848, -0.1345553 51.5102187, "
    "-0.1344994 51.510216, -0.1344519 51.51022, -0.1344024 51.5102344, "
    "-0.1343994 51.5102367, -0.1343726 51.5102578, -0.1343369 51.5103096, "
    "-0.1342521 51.510441, -0.1340647 51.5106759, -0.1339251 51.5108473, "
    "-0.1338178 51.5109648, -0.133485 51.511293, -0.1331837 51.5114941, "
    "-0.1330286 51.5116095, -0.1329784 51.5116414, -0.1328557 51.5116994, "
    "-0.1325545 51.5118259, -0.1322504 51.5119561, -0.1314392 51.5122874, "
    "-0.1304367 51.5126664, -0.1299327 51.5128501, -0.1297379 51.5129546, "
    "-0.1294894 51.5130997, -0.1293784 51.5131726, -0.1294407 51.5131871, "
    "-0.1294805 51.5132494, -0.1295057 51.5132957, -0.1295108 51.5133475, "
    "-0.1294837 51.5134495, -0.1294596 51.5134934, -0.1294072 51.5135236, "
    "-0.1294217 51.5135966, -0.1294362 51.5136317, -0.1295306 51.5137946, "
    "-0.129567 51.5138695, -0.129636 51.5139849, -0.1296451 51.5140012, "
    "-0.1296867 51.5141046, -0.1297291 51.5142121, -0.1297884 51.5143634, "
    "-0.1298403 51.5145007, -0.1299175 51.514682, -0.1300152 51.51491, "
    "-0.1301523 51.5152675, -0.1306884 51.5151383, -0.1307937 51.5153197, "
    "-0.1308895 51.5154833, -0.1307956 51.5155047, -0.1308833 51.5155937, "
    "-0.1309068 51.5155977, -0.1309316 51.5155987, -0.1309506 51.5155934, "
    "-0.1311246 51.5158882, -0.1312572 51.5160966, -0.1312752 51.5161264, "
    "-0.1309837 51.5161582, -0.1307013 51.5161832, -0.1307307 51.516365, "
    "-0.1313097 51.5163318, -0.1321982 51.5162733, -0.1329511 51.516213, "
    "-0.1335465 51.5161508, -0.1353119 51.5159597, -0.1360908 51.5158808, "
    "-0.1376583 51.5157261, -0.1386282 51.5156109, -0.1393623 51.5155356, "
    "-0.1401402 51.5154556, -0.1411432 51.5153355, -0.1415848 51.5152797, "
    "-0.1417661 51.5152281, -0.1418295 51.5150964))"
)


@pytest.fixture(scope="session")
def soho_polygon_gdf() -> gpd.GeoDataFrame:
    """Hardcoded boundary polygon for SOHO, London (no network needed)."""
    poly = wkt.loads(_SOHO_LONDON_WKT)
    return gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")


# ---------------------------------------------------------------------------
# OSM API retry helper
# ---------------------------------------------------------------------------


def call_with_osm_retry(func, *args, max_retries: int = 3, wait: int = 300, **kwargs):
    """
    Call a function and retry on RuntimeError (OSM 429/504) up to max_retries times.

    Parameters
    ----------
    func : callable
        The function to call.
    *args
        Positional arguments passed to func.
    max_retries : int
        Maximum number of attempts before re-raising.
    wait : int
        Seconds to wait between retries.
    **kwargs
        Keyword arguments passed to func.
    """
    for attempt in range(max_retries):
        try:
            return func(*args, **kwargs)
        except RuntimeError:
            if attempt == max_retries - 1:
                raise
            time.sleep(wait)
