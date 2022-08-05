import os
import numpy as np
import pandas as pd
import geopandas as gpd
import mercantile

import warnings

from shapely.ops import polygonize, unary_union
from shapely.geometry import box, Polygon, Point, LineString, MultiPolygon
from shapely.validation import make_valid

from collections import defaultdict
import h3
from sklearn.cluster import AgglomerativeClustering


def count_poi(df, points):
    """
    Counts the number POI in each tile

    Parameters
    ----------
    df : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles (polygons)
    points : geopandas.GeoDataFrame
        GeoDataFrame containing the POI

    Returns
    --------
    final_df : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles and the POI count
    """

    pointsInPolygon = gpd.sjoin(df, points, how="left", predicate="contains")
    pointsInPolygon["count"] = 1
    pointsInPolygon.reset_index(inplace=True)

    tmp_a = pointsInPolygon.groupby(by="quadkey").count()
    tmp_a = tmp_a["count"].reset_index()
    tmp_a = tmp_a.sort_values(by="quadkey", ascending=True)

    tmp_b = df.reset_index().sort_values(by="quadkey", ascending=True)

    final_df = pd.merge(tmp_a, tmp_b, on="quadkey")
    final_gdf = final_df[["quadkey", "count", "geometry", "children_id"]]
    final_gdf = gpd.GeoDataFrame(final_gdf, geometry="geometry")

    return final_gdf


##### squares


def get_squares_polyfill(gdf, zoom_level):
    """
    Square tessellation based on the quadKeys concept

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles (polygons)
    zoom_level : int
        Resolution, which controls the square sizes

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the squares
    """
    geom_name = gdf.geometry.name
    temp_dfs = []

    for idx, rows in gdf.iterrows():
        gdf_geometry = rows[geom_name]
        bbox = gdf_geometry.bounds
        tiles = mercantile.tiles(bbox[0], bbox[1], bbox[2], bbox[3], zoom_level)
        temp_rows = []
        for tile in tiles:
            temp_row = rows.copy()
            square = box(*mercantile.bounds(tile))
            if square.intersects(gdf_geometry):
                temp_row[geom_name] = square
                temp_row["quadkey"] = mercantile.quadkey(tile)

                child_ids = mercantile.children(tile)
                temp_row["children_id"] = list(
                    mercantile.quadkey(c_tile) for c_tile in child_ids
                )

                temp_rows.append(temp_row)
        temp_dfs.append(pd.DataFrame(temp_rows))

    df = pd.concat(temp_dfs)
    df = df.reset_index(drop=True)

    gdf = gpd.GeoDataFrame(df, geometry=geom_name, crs="epsg:4326")

    return gdf


##### hexagons


def get_h3_hexagons(gdf, resolution):
    """
    Hexagon tessellation based on the h3 implementation of Uber

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles (polygons)
    resolution : int
        Resolution, which controls the hexagon sizes

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the hexagons
    """
    if type(gdf.geometry.iloc[0]) == Polygon:
        hexs = h3.polyfill(
            gdf.geometry[0].__geo_interface__, resolution, geo_json_conformant=True
        )
        polygonise = lambda hex_id: Polygon(
            h3.h3_to_geo_boundary(hex_id, geo_json=True)
        )
        all_polys = gpd.GeoSeries(
            list(map(polygonise, hexs)), index=hexs, crs="EPSG:4326"
        )

        gdf = gpd.GeoDataFrame(geometry=all_polys, crs="EPSG:4326")
        return gdf

    elif type(gdf.geometry.iloc[0]) == MultiPolygon:
        parts_lst = []
        for idx, row in gdf.explode(index_parts=True).loc[0].iterrows():
            hexs = h3.polyfill(
                row.geometry.__geo_interface__, resolution, geo_json_conformant=True
            )

            polygonise = lambda hex_id: Polygon(
                h3.h3_to_geo_boundary(hex_id, geo_json=True)
            )
            all_polys = gpd.GeoSeries(
                list(map(polygonise, hexs)), index=hexs, crs="EPSG:4326"
            )

            gdf = gpd.GeoDataFrame(geometry=all_polys)

            parts_lst.append(gdf)

        return pd.concat(parts_lst)


##### adaptive squares


def get_adaptive_squares(input_gdf, threshold):
    """
    Adaptive tessellation. Subdivides all squares of input tessellation where threshold is exceeded.

    Parameters
    ----------
    input_gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the tiles (polygons)
    threshold : int
        threshold, which controls the division of squares

    Returns
    --------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing the squares
    """

    gdf = input_gdf.copy()
    # DataFrame where the squares will be subdivided
    gdf_exceeded = gdf[gdf["count"] >= threshold]

    for idx, row in gdf_exceeded.iterrows():
        children = gdf_exceeded.loc[[idx]]["children_id"].values[0]
        gdf.drop([idx], inplace=True)

        for child in children:
            new_row = row.copy()
            child_tile = mercantile.quadkey_to_tile(child)

            new_row["quadkey"] = child
            new_row["geometry"] = box(*mercantile.bounds(child_tile))
            grand_children = mercantile.children(child_tile)
            new_row["children_id"] = list(
                mercantile.quadkey(c_tile) for c_tile in grand_children
            )

            tmp_df = pd.DataFrame(new_row).transpose()
            tmp_gdf = gpd.GeoDataFrame(tmp_df, geometry="geometry", crs="epsg:4326")

            gdf = pd.concat([gdf, tmp_gdf], axis=0)

    gdf.index = gdf.reset_index(drop=True).index
    return gdf


##### voronoi diagram


def voronoi_polygons(sp_voronoi_obj, diameter):
    """
    Voronoi Diagram. Create Voronoi polygons based on scipy Voronoi object

    Parameters
    ----------
    sp_voronoi_obj : scipy.spatial.Voronoi
        scipy Voronoi object is created using the point coordinates
    diameter : float
        controls the size of infinite polygons. It should be large enough
        depending on the input values.

    Returns
    --------
    gdf : shapely.Polygon
        Voronoi polygons
    """
    centroid = sp_voronoi_obj.points.mean(axis=0)

    ridge_direction = defaultdict(list)
    for (p, q), rv in zip(sp_voronoi_obj.ridge_points, sp_voronoi_obj.ridge_vertices):
        u, v = sorted(rv)
        if u == -1:
            tangent = sp_voronoi_obj.points[q] - sp_voronoi_obj.points[p]
            normal = np.array([-tangent[1], tangent[0]]) / np.linalg.norm(tangent)
            midpoint = sp_voronoi_obj.points[[p, q]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - centroid, normal)) * normal
            ridge_direction[p, v].append(direction)
            ridge_direction[q, v].append(direction)

    for i, r in enumerate(sp_voronoi_obj.point_region):
        region = sp_voronoi_obj.regions[r]
        if -1 not in region:
            # finite regions
            yield Polygon(sp_voronoi_obj.vertices[region])
            continue

        # infinite regions
        inf_vertex_id = region.index(-1)
        previous_vertex_id = region[(inf_vertex_id - 1) % len(region)]
        next_vertex_id = region[(inf_vertex_id + 1) % len(region)]
        if previous_vertex_id == next_vertex_id:
            dir_j, dir_k = ridge_direction[i, previous_vertex_id]
        else:
            (dir_j,) = ridge_direction[i, previous_vertex_id]
            (dir_k,) = ridge_direction[i, next_vertex_id]

        length = 2 * diameter / np.linalg.norm(dir_j + dir_k)

        finite_part = sp_voronoi_obj.vertices[
            region[inf_vertex_id + 1 :] + region[:inf_vertex_id]
        ]
        extra_edge = [
            sp_voronoi_obj.vertices[previous_vertex_id] + dir_j * length,
            sp_voronoi_obj.vertices[next_vertex_id] + dir_k * length,
        ]
        yield Polygon(np.concatenate((finite_part, extra_edge)))


##### city blocks


def split_linestring(df):
    """
    The linestring, which are split up will have the same osmid
    since it's still the same street

    Parameters
    ----------
    df: pandas.DataFrame
        dataframe of shapely.LineStrings containing more than two points

    Returns
    --------
    dataset: pandas.DataFrame
        dataframe with shapely.linestring with two points each
    """
    linestrings = []
    osmid = []

    for idx, row in df.iterrows():
        if len(row["geometry"].coords) == 2:
            linestrings.append(row["geometry"])
            osmid.append(row["osmid"])
        else:
            for i in range(0, len(row["geometry"].coords) - 1):
                p1 = Point(row["geometry"].coords[i][0], row["geometry"].coords[i][1])
                p2 = Point(
                    row["geometry"].coords[i + 1][0], row["geometry"].coords[i + 1][1]
                )
                linestrings.append(LineString([p1, p2]))
                osmid.append(row["osmid"])

    dataset = gpd.GeoDataFrame({"osmid": osmid, "geometry": linestrings})
    return dataset


def explode(gdf):
    """
    Multipolygon into multiple Polygons

    Parameters
    ----------
    gdf:  geopandas.GeoDataFrame
        GeoDataFrame that can have multi-part geometries that will be exploded

    Returns
    --------
    gdf_out: geopandas.GeoDataFrame
        GeoDataFrame with single geometries
    """
    gs = gdf.explode(index_parts=True)
    gdf2 = gs.reset_index().rename(columns={0: "geometry"})
    gdf_out = gdf2.merge(
        gdf.drop("geometry", axis=1),
        left_on="level_0",
        right_index=True,
    )
    gdf_out = gdf_out.set_index(["level_0", "level_1"]).set_geometry("geometry")
    gdf_out.crs = gdf.crs
    return gdf_out


def get_hierarchical_clustering_parameter(coordinates, threshold):
    """
    This function returns the distance_threshold that is used in
    the hierarchical clustering algorithm. Therefore, different
    distant_threshold are tested until the hierarchical clustering
    return less clusters than the input threshold. If no dist_threshold
    fulfills the threshold-requirement nothing is returned. In the
    main method tesspy.Tessellation.city_blocks() is a differentiation
    between both cases. If th is an integer, the hierarchical clustering
    use distance_threshold parameter, if th is None, the hierarchical
    clustering use n_clusters parameter with the users input
    n_polygons

    Parameters
    ----------
    coordinates: numpy.ndarray
        stock of coordinates of Data
    threshold: int
        positive number indicating the number of Polygons that should not be exceeded

    Returns
    --------
    th: int
        defines the distance_threshold for hierarchical clustering
    """
    dist_threshold = [i * 100 for i in range(2, 13)]
    for th in dist_threshold:
        model = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=th,
            affinity="euclidean",
            compute_full_tree=True,
        )
        model.fit(coordinates)
        labels = model.labels_
        nb_clusters = len(set(labels)) - (1 if -1 in labels else 0)

        if nb_clusters < threshold:
            return th


def create_blocks(road_network):
    """
    This function uses polygonize to create blocks by using road data.

    Parameters
    ----------
    road_network: geopandas.GeoDataFrame
        GeoDataFrame (coming from RoadData class) containing street segments

    Returns
    --------
    blocks: geopandas.GeoDataFrame
        GeoDataFrame with polygons, that were created by using the road data
    """
    if hasattr(road_network, "geometry"):
        block_faces = list(polygonize(road_network["geometry"]))
        blocks = gpd.GeoDataFrame(geometry=block_faces).set_crs("EPSG:4326")
        return blocks
    else:
        raise AttributeError("road network data need geometry attribute!")


def get_rest_polygon(blocks, area):
    # TODO: think about making this function adapt better to islands
    #  groups with ocean and multiple polygons because now the convex
    #  hull is used. That results in ocean parts becoming city blocks
    """
    This function creates the rest polygons that can occur by creating
    city blocks. Dead ends or vegetation at the boundary of the area will
    not be defined as blocks because there are no road to define a block.
    Since the tessellation should cover the whole area those "rest polygons"
    will be created by subtracting all blocks from the area to fill these gaps.
    The rest_polygons are most likely a Multi-Polygon which is exploded into many polygons.
    Using this, the city blocks method tessellates the whole area without gaps.

    Parameters
    ----------
    blocks: geopandas.GeoDataFrame
        GeoDataFrame containing city blocks
    area: geopandas.GeoDataFrame
        boundary polygon

    Returns
    --------
    rest_polygons: geopandas.GeoDataFrame
        GeoDataFrame containing rest polygons
    """

    if hasattr(blocks, "geometry") and hasattr(area, "geometry"):

        blocks["geometry"] = blocks["geometry"].apply(lambda x: make_valid(x))

        merged_polygons = gpd.GeoSeries(unary_union(blocks["geometry"].values))
        merged_polygons.set_crs("EPSG:4326", allow_override=True, inplace=True)

        rest = area.difference(merged_polygons)
        rest = gpd.GeoDataFrame(rest)
        rest = rest.rename(columns={0: "geometry"}).set_geometry("geometry")

        rest_polygons = explode(rest)
        rest_polygons.reset_index(inplace=True)
        rest_polygons.drop(columns=["level_0"], inplace=True)
        rest_polygons.rename(columns={"level_1": "osm_id"}, inplace=True)

        return rest_polygons

    else:
        raise ValueError(
            "Initial defined city blocks and the area need a geometry attribute."
        )
