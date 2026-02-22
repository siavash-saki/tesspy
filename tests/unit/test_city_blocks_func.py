"""
Unit tests for city block functions (no network required).

Tests for the new ``create_city_blocks`` / ``merge_city_blocks`` API and
deprecation-warning tests for the legacy helpers.
"""

import warnings

import geopandas as gpd
import pytest
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from tesspy.methods.city_blocks import (
    create_blocks,
    create_city_blocks,
    explode,
    get_rest_polygon,
    merge_city_blocks,
    split_linestring,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def square_area():
    """A 3x3 square study area."""
    poly = Polygon([(0, 0), (3, 0), (3, 3), (0, 3)])
    return gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")


@pytest.fixture()
def grid_roads():
    """A grid of roads that divides a 3x3 square into 6 blocks (3 cols x 2 rows)."""
    lines = [
        LineString([(0, 0), (3, 0)]),
        LineString([(0, 0), (0, 3)]),
        LineString([(3, 0), (3, 3)]),
        LineString([(0, 3), (3, 3)]),
        LineString([(1, 0), (1, 3)]),
        LineString([(2, 0), (2, 3)]),
        LineString([(0, 1.5), (3, 1.5)]),
    ]
    return gpd.GeoDataFrame(geometry=lines, crs="EPSG:4326")


# ---------------------------------------------------------------------------
# create_city_blocks tests
# ---------------------------------------------------------------------------


class TestCreateCityBlocks:
    def test_basic_block_count(self, grid_roads, square_area):
        """Grid roads inside a square produce 6 blocks (3 cols x 2 rows)."""
        blocks = create_city_blocks(grid_roads, square_area)
        assert isinstance(blocks, gpd.GeoDataFrame)
        assert len(blocks) == 6

    def test_full_coverage(self, grid_roads, square_area):
        """Union of all blocks should equal the study area (no gaps)."""
        blocks = create_city_blocks(grid_roads, square_area)
        block_union = unary_union(blocks.geometry)
        area_polygon = square_area.geometry.iloc[0]
        sym_diff_area = block_union.symmetric_difference(area_polygon).area
        assert sym_diff_area < 1e-10

    def test_all_polygons(self, grid_roads, square_area):
        """All output geometries should be single Polygons."""
        blocks = create_city_blocks(grid_roads, square_area)
        assert (blocks.geom_type == "Polygon").all()

    def test_blocks_inside_area(self, grid_roads, square_area):
        """Every block's representative point should be inside the area."""
        blocks = create_city_blocks(grid_roads, square_area)
        area_polygon = square_area.geometry.iloc[0]
        inside = blocks.geometry.representative_point().within(area_polygon)
        assert inside.all()

    def test_crs_set(self, grid_roads, square_area):
        blocks = create_city_blocks(grid_roads, square_area)
        assert blocks.crs is not None

    def test_empty_road_network(self, square_area):
        """With no roads, the whole area becomes one block."""
        empty_roads = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
        blocks = create_city_blocks(empty_roads, square_area)
        assert len(blocks) == 1

    def test_roads_extending_beyond_area(self, square_area):
        """Roads extending beyond the area should not produce outside blocks."""
        lines = [
            LineString([(-1, 1.5), (4, 1.5)]),
            LineString([(1.5, -1), (1.5, 4)]),
        ]
        roads = gpd.GeoDataFrame(geometry=lines, crs="EPSG:4326")
        blocks = create_city_blocks(roads, square_area)
        area_polygon = square_area.geometry.iloc[0]
        inside = blocks.geometry.representative_point().within(area_polygon)
        assert inside.all()
        assert len(blocks) == 4

    def test_with_real_fixture_data(self, road_data_lille, soho_polygon_gdf):
        """Smoke test with pre-downloaded Lille road data."""
        blocks = create_city_blocks(road_data_lille, soho_polygon_gdf)
        assert isinstance(blocks, gpd.GeoDataFrame)
        assert len(blocks) > 0
        assert blocks.crs is not None


# ---------------------------------------------------------------------------
# merge_city_blocks tests
# ---------------------------------------------------------------------------


class TestMergeCityBlocks:
    def test_reduces_count(self, grid_roads, square_area):
        """Merging should reduce block count to the target."""
        blocks = create_city_blocks(grid_roads, square_area)
        merged = merge_city_blocks(blocks, n_polygons=3)
        assert len(merged) == 3

    def test_preserves_coverage(self, grid_roads, square_area):
        """Merged blocks should cover the same total area."""
        blocks = create_city_blocks(grid_roads, square_area)
        original_area = unary_union(blocks.geometry).area
        merged = merge_city_blocks(blocks, n_polygons=2)
        merged_area = unary_union(merged.geometry).area
        assert abs(original_area - merged_area) < 1e-10

    def test_contiguous_results(self, grid_roads, square_area):
        """With adjacency constraint, merged groups should be contiguous Polygons."""
        blocks = create_city_blocks(grid_roads, square_area)
        merged = merge_city_blocks(blocks, n_polygons=3)
        assert (merged.geom_type == "Polygon").all()


# ---------------------------------------------------------------------------
# Deprecation warning tests for legacy functions
# ---------------------------------------------------------------------------


class TestDeprecatedFunctions:
    def test_split_linestring_warns(self, road_data_lille):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            split_linestring(road_data_lille)
        dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert len(dep) >= 1
        assert "split_linestring" in str(dep[0].message)

    def test_split_linestring_still_works(self, road_data_lille):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = split_linestring(road_data_lille)
        assert isinstance(result, gpd.GeoDataFrame)
        assert len(result) > len(road_data_lille)

    def test_explode_warns(self):
        gdf = gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])], crs="EPSG:4326"
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            explode(gdf)
        dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert len(dep) >= 1
        assert "explode" in str(dep[0].message)

    def test_create_blocks_warns(self, road_data_lille):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            create_blocks(road_data_lille)
        dep = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert len(dep) >= 1
        assert "create_blocks" in str(dep[0].message)

    def test_get_rest_polygon_warns(self, soho_polygon_gdf):
        blocks_gdf = gpd.GeoDataFrame(
            geometry=[Polygon([(0, 0), (0.001, 0), (0.001, 0.001), (0, 0.001)])],
            crs="EPSG:4326",
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            get_rest_polygon(blocks_gdf, soho_polygon_gdf)
        dep = [
            w
            for w in caught
            if issubclass(w.category, DeprecationWarning)
            and "get_rest_polygon" in str(w.message)
        ]
        assert len(dep) >= 1
