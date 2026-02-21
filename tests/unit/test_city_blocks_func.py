"""
Unit tests for city block helper functions (no network required).
Uses pre-downloaded local fixture data (road_data_Lille.pkl).
"""

import geopandas as gpd
import pytest

from tesspy.methods.city_blocks import (
    create_blocks,
    explode,
    get_rest_polygon,
    split_linestring,
)


def test_split_linestring(road_data_lille):
    """Test that split_linestring produces more 2-point segments than original."""
    split_data = split_linestring(road_data_lille)

    assert isinstance(split_data, gpd.GeoDataFrame)
    assert len(split_data) > len(road_data_lille)
    assert "geometry" in split_data.columns
    assert "osmid" in split_data.columns

    # All resulting geometries should have exactly 2 coordinate points
    coord_counts = split_data["geometry"].apply(lambda g: len(list(g.coords)))
    assert (coord_counts == 2).all()


def test_create_blocks(road_data_lille):
    """Test that create_blocks produces closed polygon blocks from road data."""
    split_data = split_linestring(road_data_lille)
    blocks = create_blocks(split_data)

    assert isinstance(blocks, gpd.GeoDataFrame)
    assert len(blocks) > 0
    assert blocks.crs is not None


def test_explode():
    """Test that explode splits MultiPolygons into individual Polygons."""
    from shapely.geometry import MultiPolygon, Polygon

    poly1 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    poly2 = Polygon([(2, 0), (3, 0), (3, 1), (2, 1)])
    multi = MultiPolygon([poly1, poly2])

    gdf = gpd.GeoDataFrame(geometry=[multi], crs="EPSG:4326")
    result = explode(gdf)

    assert isinstance(result, gpd.GeoDataFrame)
    assert len(result) == 2
