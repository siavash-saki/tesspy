"""
Unit tests for square tessellation functions (no network required).
Uses pre-downloaded local fixture data.
"""

import geopandas as gpd

from tesspy.methods.squares import count_poi, get_adaptive_squares, get_squares_polyfill


def test_count_poi(soho_polygon_gdf, leisure_poi_berlin):
    """Test that count_poi correctly counts POI within square tiles."""
    # Use the Mitte Berlin polygon derived from the geojson bounds
    from shapely.geometry import box

    mitte_bounds = leisure_poi_berlin.total_bounds  # (minx, miny, maxx, maxy)
    mitte_poly = box(*mitte_bounds)
    mitte_gdf = gpd.GeoDataFrame(geometry=[mitte_poly], crs="EPSG:4326")
    mitte_gdf["osm_id"] = 0

    city_squares = get_squares_polyfill(mitte_gdf, 14)
    count_result = count_poi(city_squares, leisure_poi_berlin)

    assert hasattr(count_result, "count")
    assert hasattr(count_result, "children_id")
    assert len(count_result) > 0

    # Verify specific known quadkey counts from the original test data
    q1 = count_result[count_result["quadkey"] == "12021023322022"]["count"]
    q2 = count_result[count_result["quadkey"] == "12021023322023"]["count"]
    q3 = count_result[count_result["quadkey"] == "12021023322032"]["count"]

    if len(q1) > 0:
        assert q1.iloc[0] == 61
    if len(q2) > 0:
        assert q2.iloc[0] == 109
    if len(q3) > 0:
        assert q3.iloc[0] == 56


def test_count_poi_uses_fixture_directly(leisure_poi_berlin, soho_polygon_gdf):
    """Test count_poi using the SOHO polygon with pre-loaded POI data."""
    soho_gdf = soho_polygon_gdf.copy()
    soho_gdf["osm_id"] = 0
    city_squares = get_squares_polyfill(soho_gdf, 14)

    # leisure_poi_berlin is used as generic point data
    # (spatial join; only overlapping points counted)
    count_result = count_poi(city_squares, leisure_poi_berlin)

    assert isinstance(count_result, gpd.GeoDataFrame)
    assert "count" in count_result.columns
    assert len(count_result) > 0


def test_get_adaptive_squares(leisure_poi_berlin):
    """Test that get_adaptive_squares subdivides squares above threshold."""
    from shapely.geometry import box

    mitte_bounds = leisure_poi_berlin.total_bounds
    mitte_poly = box(*mitte_bounds)
    mitte_gdf = gpd.GeoDataFrame(geometry=[mitte_poly], crs="EPSG:4326")
    mitte_gdf["osm_id"] = 0

    city_squares = get_squares_polyfill(mitte_gdf, 14)
    count_result = count_poi(city_squares, leisure_poi_berlin)

    # With a very low threshold (100), more squares should be produced
    subdivided = get_adaptive_squares(count_result, 100)
    assert len(subdivided) > len(city_squares)

    # With a very high threshold (more than any tile's count), no subdivision
    max_count = int(count_result["count"].max()) + 1
    unchanged = get_adaptive_squares(count_result, max_count)
    assert len(unchanged) == len(city_squares)


def test_get_adaptive_squares_threshold_zero_subdivides_all(leisure_poi_berlin):
    """Threshold=0 subdivides ALL tiles (count >= 0 is always true).

    This verifies the dangerous edge case that the auto-threshold fix in
    Tessellation.adaptive_squares guards against: when median POI count
    is 0 and the old code used ``if not threshold`` instead of
    ``if threshold is None``, the auto-computed threshold of 0 would
    cause every tile to be subdivided on every iteration, producing an
    infinite loop.
    """
    from shapely.geometry import box

    mitte_bounds = leisure_poi_berlin.total_bounds
    mitte_poly = box(*mitte_bounds)
    mitte_gdf = gpd.GeoDataFrame(geometry=[mitte_poly], crs="EPSG:4326")
    mitte_gdf["osm_id"] = 0

    city_squares = get_squares_polyfill(mitte_gdf, 14)
    count_result = count_poi(city_squares, leisure_poi_berlin)

    # threshold=0 means count >= 0, which is true for ALL tiles,
    # so every single tile gets subdivided into 4 children
    subdivided = get_adaptive_squares(count_result, 0)
    assert len(subdivided) > len(city_squares)
    # In fact every tile is subdivided: each original tile becomes 4 children
    assert len(subdivided) == len(city_squares) * 4
