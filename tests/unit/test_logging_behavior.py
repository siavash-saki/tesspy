"""
Unit tests for runtime logging behavior across tesspy modules.
"""

import logging

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, Polygon

import tesspy.data._geo as geo_module
import tesspy.data.poi as poi_module
import tesspy.data.roads as roads_module
import tesspy.tessellation as tess_module
from tesspy.data.poi import POIdata
from tesspy.data.roads import RoadData
from tesspy.tessellation import Tessellation


def _fake_features_from_polygon(*args, **kwargs):
    """Return a minimal GeoDataFrame mimicking osmnx.features_from_polygon."""
    return gpd.GeoDataFrame(
        {"amenity": ["cafe", "school"]},
        geometry=[Point(-0.13, 51.51), Point(-0.131, 51.512)],
        crs="EPSG:4326",
    )


def test_poi_verbose_true_emits_progress_events(soho_polygon_gdf, monkeypatch, caplog):
    monkeypatch.setattr(
        poi_module.ox, "features_from_polygon", _fake_features_from_polygon
    )

    poi = POIdata(
        soho_polygon_gdf,
        poi_categories=["amenity"],
        timeout=60,
        verbose=True,
    )
    with caplog.at_level(logging.DEBUG, logger="tesspy.data.poi"):
        result = poi.get_poi_data()

    messages = [record.getMessage() for record in caplog.records]
    assert len(result) > 0
    assert any("event=poi.fetch.start" in msg for msg in messages)
    assert any("event=poi.parse.done" in msg for msg in messages)


def test_poi_verbose_false_suppresses_progress_events(
    soho_polygon_gdf, monkeypatch, caplog
):
    monkeypatch.setattr(
        poi_module.ox, "features_from_polygon", _fake_features_from_polygon
    )

    poi = POIdata(
        soho_polygon_gdf,
        poi_categories=["amenity"],
        timeout=60,
        verbose=False,
    )
    with caplog.at_level(logging.DEBUG, logger="tesspy.data.poi"):
        poi.get_poi_data()

    messages = [record.getMessage() for record in caplog.records]
    assert not any("event=poi.fetch.start" in msg for msg in messages)
    assert not any("event=poi.parse.done" in msg for msg in messages)


def test_poi_empty_response_logs_warning(soho_polygon_gdf, monkeypatch, caplog):
    empty_gdf = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    monkeypatch.setattr(
        poi_module.ox,
        "features_from_polygon",
        lambda *a, **kw: empty_gdf,
    )

    poi = POIdata(
        soho_polygon_gdf,
        poi_categories=["amenity"],
        timeout=60,
        verbose=False,
    )

    with caplog.at_level(logging.DEBUG, logger="tesspy.data.poi"):
        with pytest.raises(ValueError):
            poi.get_poi_data()

    log_records = [
        r for r in caplog.records if "event=poi.fetch.empty" in r.getMessage()
    ]
    assert len(log_records) == 1
    assert log_records[0].levelno == logging.WARNING


def test_roaddata_emits_progress_events(soho_polygon_gdf, monkeypatch, caplog):
    class DummyGraph:
        def to_undirected(self):
            return self

    monkeypatch.setattr(
        roads_module.ox,
        "graph_from_polygon",
        lambda *args, **kwargs: DummyGraph(),
    )
    monkeypatch.setattr(roads_module.ox, "project_graph", lambda graph, to_crs: graph)
    monkeypatch.setattr(
        roads_module.ox,
        "graph_to_gdfs",
        lambda *args, **kwargs: gpd.GeoDataFrame(
            {"osmid": [1], "highway": ["residential"]},
            geometry=[LineString([(0, 0), (1, 1)])],
            crs="EPSG:4326",
        ),
    )

    with caplog.at_level(logging.DEBUG, logger="tesspy.data.roads"):
        RoadData(soho_polygon_gdf, detail_deg=3, verbose=True).get_road_network()

    messages = [record.getMessage() for record in caplog.records]
    assert any("event=roads.filter.selected" in msg for msg in messages)
    assert any("event=roads.fetch.start" in msg for msg in messages)
    assert any("event=roads.fetch.done" in msg for msg in messages)


def test_tessellation_logs_cache_hit_events(soho_polygon_gdf, monkeypatch, caplog):
    tess = Tessellation(soho_polygon_gdf)
    tess.poi_dataframe = pd.DataFrame(
        {
            "center_longitude": [-0.13],
            "center_latitude": [51.51],
            "amenity": [True],
        }
    )
    tess.queried_highway_types = tess.osm_highway_types()
    tess.road_network = gpd.GeoDataFrame(
        {"osmid": [1]},
        geometry=[LineString([(0, 0), (1, 1)])],
        crs="EPSG:4326",
    )

    area_poly = soho_polygon_gdf.geometry.iloc[0]

    monkeypatch.setattr(
        tess_module,
        "get_squares_polyfill",
        lambda area, res: gpd.GeoDataFrame(
            {"osm_id": [1], "children_id": [[]], "quadkey": ["q1"]},
            geometry=[area.geometry.iloc[0]],
            crs="EPSG:4326",
        ),
    )
    monkeypatch.setattr(
        tess_module,
        "count_poi",
        lambda df, points: df.assign(count=1),
    )
    monkeypatch.setattr(
        tess_module.gpd,
        "sjoin",
        lambda left, right, *args, **kwargs: left.assign(index_right=0),
    )
    monkeypatch.setattr(
        tess_module,
        "create_city_blocks",
        lambda roads, area: gpd.GeoDataFrame(geometry=[area_poly], crs="EPSG:4326"),
    )
    monkeypatch.setattr(tess_module, "_check_valid_geometry_gdf", lambda gdf: gdf)

    with caplog.at_level(logging.DEBUG, logger="tesspy.tessellation"):
        tess.adaptive_squares(
            start_resolution=10,
            poi_categories=["amenity"],
            threshold=1,
            verbose=True,
        )
        tess.city_blocks(
            n_polygons=None,
            detail_deg=None,
            verbose=True,
        )

    messages = [record.getMessage() for record in caplog.records]
    assert any("event=poi.cache.hit" in msg for msg in messages)
    assert any("event=roads.cache.hit" in msg for msg in messages)


def test_count_poi_per_tile_forwards_verbose(soho_polygon_gdf, monkeypatch):
    captured = {"verbose": None}

    class FakePOIdata:
        def __init__(self, area, poi_categories, timeout, verbose):
            captured["verbose"] = verbose
            self.poi_categories = poi_categories

        def get_poi_data(self):
            return pd.DataFrame(
                {
                    "center_longitude": [0.5],
                    "center_latitude": [0.5],
                    "amenity": [True],
                }
            )

    monkeypatch.setattr(geo_module, "POIdata", FakePOIdata)
    monkeypatch.setattr(
        geo_module.gpd,
        "sjoin",
        lambda left, right, *args, **kwargs: pd.DataFrame(
            {"quadkey": [left["quadkey"].iloc[0]], "value": ["amenity"]}
        ),
    )

    tile = gpd.GeoDataFrame(
        {"quadkey": ["q1"]},
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        crs="EPSG:4326",
    )
    result = geo_module.count_poi_per_tile(
        soho_polygon_gdf,
        tile,
        poi_categories=["amenity"],
        verbose=True,
    )

    assert captured["verbose"] is True
    assert "amenity" in result.columns
