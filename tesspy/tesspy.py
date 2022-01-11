import hdbscan
import numpy as np
import pandas as pd
import geopandas as gpd
from geopandas import sjoin
import osmnx as ox
import h3
from babelgrid import Babel
import shapely
from shapely.geometry import Point, Polygon, LineString, mapping, MultiPoint
from shapely.ops import polygonize, cascaded_union
from sklearn.cluster import AgglomerativeClustering, KMeans
from scipy.spatial import Voronoi
from collections import defaultdict


def get_admin_polygon(city: str):
    """
    :param city:
    :return: simple GeoDataFrame containing only the Polygon of the city
    """
    try:
        df_city = ox.geocode_to_gdf(city)
        df_city = df_city[["osm_id", "geometry"]]
        df_city = df_city.rename({"osm_id": "osmid"})
        return df_city

    except:
        print("Input must be a city in string format")


def get_bbox(city: str):
    """

    :param city:
    :return: Array of Int which represent the bounding boy for the given city in the format
    [north, south,east,west]
    """
    gdf_city = ox.geocode_to_gdf(city)
    bbox = [gdf_city["bbox_north"].values[0],
            gdf_city["bbox_south"].values[0],
            gdf_city["bbox_east"].values[0],
            gdf_city["bbox_west"].values[0]]
    return bbox


def split_linestring(df):
    """

    :param df:
    :return: DataFrame with shapely.linestring with two points each
    The linestring, which are split up will have the same osmid since its still the same street
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
                p2 = Point(row["geometry"].coords[i + 1][0], row["geometry"].coords[i + 1][1])
                linestrings.append(LineString([p1, p2]))
                osmid.append(row["osmid"])

    dataset = gpd.GeoDataFrame({"osmid": osmid, "geometry": linestrings})
    return dataset


def explode(gdf):
    """

    :param gdf: which can have multi-part geometries that will be exploded
    :return: GeoDataFrame with single geometries
    example: Multipolygon -> multiple Polygons
    """
    gs = gdf.explode(index_parts=True)
    gdf2 = gs.reset_index().rename(columns={0: "geometry"})
    gdf_out = gdf2.merge(gdf.drop("geometry", axis=1), left_on="level_0", right_index=True, )
    gdf_out = gdf_out.set_index(["level_0", "level_1"]).set_geometry("geometry")
    gdf_out.crs = gdf.crs
    return gdf_out


def count_poi(df, points):
    """
    :param df: Geodataframe of LGUs
    :param points: GeoDataFrame of POI
    :return: GeoDataFrame of LGUs with "count" column that counts points per LGU
    """

    pointsInPolygon = gpd.sjoin(df, points, how="left", predicate='contains')
    pointsInPolygon['count'] = 1
    pointsInPolygon.reset_index(inplace=True)
    tmp_a = pointsInPolygon.groupby(by='tile_id').count()['count'].reset_index().sort_values(by="tile_id",
                                                                                             ascending=True)
    tmp_b = df.reset_index().sort_values(by="tile_id", ascending=True)
    final_df = pd.merge(tmp_a, tmp_b, on="tile_id")
    final_df.set_index("tile_id", inplace=True)

    return final_df


def get_h3(df, resolution):
    """
    :param df: should the geodatafram of the boundary polygon
    :param resolution: resolution for uber's h3 hexagon grid
    :return: GeoDataFrame with tessellation
    """
    area_json = mapping(df.geometry.iloc[0])
    h3_index = h3.polyfill(area_json, resolution)
    h3_polygons = [Polygon(h3.h3_to_geo_boundary(h3_idx)) for h3_idx in h3_index]
    gdf = gpd.GeoDataFrame(geometry=h3_polygons, crs='EPSG:4326')
    return gdf


def get_LGU_threshold(city_admin_boundary, res):
    """
    :param city_admin_boundary: boundary polygon of the city
    :param res: resolution for h3 tessellation
    :return: number of LGU
    """
    hexagons = city_admin_boundary.h3.polyfill_resample(res)
    return hexagons.shape[0]


def get_hdscan_parameter(coordinates, threshold):
    dist_threshold = [i * 50 for i in range(2, 13)]
    for th in dist_threshold:
        model = AgglomerativeClustering(n_clusters=None, distance_threshold=th, affinity="euclidean",
                                        compute_full_tree=True)
        model.fit(coordinates)
        labels = model.labels_
        nb_clusters = len(set(labels)) - (1 if -1 in labels else 0)

        if nb_clusters < threshold:
            return th


def adaptive_tessellation(gdf, threshold):
    """
    :param gdf: GeoDataFrame from babel (with children_id column)
    :param threshold: Threshold of max data points per LGU allowed
    :return: divided GeoDataFrame where LGUs are subdivided that exceed the threshold
    """

    # DataFrame where the squares will be subdivided
    gdf_exceeded = gdf[gdf["count"] >= threshold]

    for idx in gdf_exceeded.index:
        # get all children (so array of length 4)
        children = gdf_exceeded.loc[[str(idx)]]["children_id"].values[0]

        # Delete the row in the overall dataframe
        gdf.drop([idx], inplace=True)

        for child in children:
            child_dict = Babel('bing').id_to_tile(child).to_dict()
            child_gdf = gpd.GeoDataFrame.from_dict(child_dict, orient='index')
            child_gdf = child_gdf.transpose()
            child_gdf = child_gdf.set_geometry("shapely")
            child_gdf = child_gdf.set_crs("EPSG:4326")
            child_gdf.set_index('tile_id', inplace=True)
            gdf = gdf.append(child_gdf)

    return gdf


def voronoi_polygons(voronoi, diameter):
    """Generate shapely.geometry.Polygon objects corresponding to the
    regions of a scipy.spatial.Voronoi object, in the order of the
    input points. The polygons for the infinite regions are large
    enough that all points within a distance 'diameter' of a Voronoi
    vertex are contained in one of the infinite polygons.

    """
    centroid = voronoi.points.mean(axis=0)

    # Mapping from (input point index, Voronoi point index) to list of
    # unit vectors in the directions of the infinite ridges starting
    # at the Voronoi point and neighbouring the input point.
    ridge_direction = defaultdict(list)
    for (p, q), rv in zip(voronoi.ridge_points, voronoi.ridge_vertices):
        u, v = sorted(rv)
        if u == -1:
            # Infinite ridge starting at ridge point with index v,
            # equidistant from input points with indexes p and q.
            t = voronoi.points[q] - voronoi.points[p]  # tangent
            n = np.array([-t[1], t[0]]) / np.linalg.norm(t)  # normal
            midpoint = voronoi.points[[p, q]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - centroid, n)) * n
            ridge_direction[p, v].append(direction)
            ridge_direction[q, v].append(direction)

    for i, r in enumerate(voronoi.point_region):
        region = voronoi.regions[r]
        if -1 not in region:
            # Finite region.
            yield Polygon(voronoi.vertices[region])
            continue
        # Infinite region.
        inf = region.index(-1)  # Index of vertex at infinity.
        j = region[(inf - 1) % len(region)]  # Index of previous vertex.
        k = region[(inf + 1) % len(region)]  # Index of next vertex.
        if j == k:
            # Region has one Voronoi vertex with two ridges.
            dir_j, dir_k = ridge_direction[i, j]
        else:
            # Region has two Voronoi vertices, each with one ridge.
            dir_j, = ridge_direction[i, j]
            dir_k, = ridge_direction[i, k]

        # Length of ridges needed for the extra edge to lie at least
        # 'diameter' away from all Voronoi vertices.
        length = 2 * diameter / np.linalg.norm(dir_j + dir_k)

        # Polygon consists of finite part plus an extra edge.
        finite_part = voronoi.vertices[region[inf + 1:] + region[:inf]]
        extra_edge = [voronoi.vertices[j] + dir_j * length,
                      voronoi.vertices[k] + dir_k * length]
        yield Polygon(np.concatenate((finite_part, extra_edge)))


class TessObj:

    def __init__(self, city, resolution):
        self.resolution = resolution
        self.city = city
        pass

    def quadKey(self, city, resolution: int):
        """
        :param resolution: resolution for mircosoft bing. Int ∈ [1,...,23]
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :return:GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")
        tiles = Babel('bing').polyfill(df_city.geometry.unary_union, resolution=resolution)
        df_qk = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely', crs="EPSG:4326")
        return df_qk

    def hexagon(self, city, resolution: int):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param resolution: resolution for mircosoft bing. Int ∈ [0,...,15]
        :return: GeoDataFrame with hexagons for the given city, geometry in shapely-Polygon
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")

        df_h3 = get_h3(df_city, resolution)
        return df_h3

    def adaptive_quadkey(self, city, poi_data, start_resolution):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param poi_data: GeoDataFrame of Data that will be used as subdivide criteria
        :param start_resolution: start resolution from which the algorithm will subdivide the area
        :return: GeoDataFrame of the LGUs with additional column "count" that counts the poi_data per LGU
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")

        tiles = Babel("bing").polyfill(df_city.geometry.unary_union, resolution=start_resolution)
        df_aqk = gpd.GeoDataFrame([t.to_dict() for t in tiles], geometry='shapely')
        df_aqk = df_aqk.set_crs("EPSG:4326")
        df_aqk.set_index("tile_id", inplace=True)

        aqk_count = count_poi(df_aqk, poi_data)
        aqk_count = gpd.GeoDataFrame(aqk_count, geometry="shapely")
        threshold = int(np.median(aqk_count["count"].values))
        while max(aqk_count["count"].values) > threshold:
            df_temp = adaptive_tessellation(aqk_count, threshold)
            df_temp.drop(columns=["count"], inplace=True)

            df_temp2 = count_poi(df_temp, poi_data)
            aqk_count = gpd.GeoDataFrame(df_temp2, geometry="shapely")

        return aqk_count

    def voronoi(self, city, data, cluster_algo="k-means"):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :param data: GeoDataFrame with a geometry column representing data like POI or generators for Voronoi Diagrams
        :param cluster_algo: if not generators are used directly the spatial data is going to be clustered. specify the
        cluster algorithm
        :return: GeoDataFrame with Polygons representing the voronoi polygons
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            pass
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")

        data_locs = np.column_stack([data.geometry.x, data.geometry.y])

        # TODO: maybe include more cluster algorithms like DBSCAN,..
        if cluster_algo == "k-means":
            # TODO: Think about automating nb of clusters
            clusterer = KMeans(n_clusters=500).fit(data_locs)
            data["ClusterIDs"] = clusterer.labels_

        elif cluster_algo == "hdbscan":
            # TODO: Think about automating min_cluster_size e.g. x% of the initial data shape
            clusterer = hdbscan.HDBSCAN(min_cluster_size=15, prediction_data=True).fit(data_locs)
            soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
            data["ClusterIDs"] = [np.argmax(x) for x in soft_clusters]

        elif cluster_algo == "None":
            data["ClusterIDs"] = data.geometry
        else:
            raise ValueError("Please use one of the implemented clustering: k-mean, hdbscan, dbscan")

        convex_hulls = [MultiPoint(data[data["ClusterIDs"] == ID]["geometry"].values).convex_hull for ID in
                        data["ClusterIDs"].unique()]
        gdf_generators = gpd.GeoDataFrame(geometry=convex_hulls, crs="EPSG:4326")
        generators = gdf_generators["geometry"].centroid

        voronoi_dia = Voronoi(np.column_stack([generators.x, generators.y]))
        voronoi_poly = gpd.GeoDataFrame(geometry=[p for p in voronoi_polygons(voronoi_dia, 0.1)], crs="EPSG:4326")
        vor_polygons = voronoi_poly.intersection(df_city.geometry.iloc[0])

        return gpd.GeoDataFrame(geometry=vor_polygons)

    def cityblocks(self, city):
        """
        :param city: Must be a shapely.Polygon, GeoDataFrame Polygon or String e.g. "Frankfurt am Main"
        :return: GeoDataFrame with all cityblocks_tiles
        """
        if type(city) == shapely.geometry.polygon.Polygon:
            df_city = gpd.GeoDataFrame(geometry=[city])
            box = city.geometry.values[0].minimum_rotated_rectangle
            x, y = box.exterior.coords.xy
            bbox = [max(y), min(y), max(x), min(x)]
        elif type(city) == str:
            df_city = get_admin_polygon(city)
            bbox = get_bbox(city)
            pass
        elif type(city) == gpd.GeoDataFrame:
            df_city = city.copy()
            try:
                [df_city["bbox_north"].values[0], df_city["bbox_south"].values[0], df_city["bbox_east"].values[0],
                 df_city["bbox_west"].values[0]]
            except:
                box = city.geometry.values[0].minimum_rotated_rectangle
                x, y = box.exterior.coords.xy
                bbox = [max(y), min(y), max(x), min(x)]
            pass
        else:
            raise TypeError("City must be in format: Shapely Polygon, GeoDataFrame or String")

        cf = "['highway'~'motorway|trunk|motorway_link|trunk|trunk_link|primary|secondary|" \
             "tertiary|residential|footway|living_street|pedestrian|track|footway|path|service|cycleway']"
        G = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3], custom_filter=cf)
        G_projected = ox.project_graph(G, to_crs='epsg:4326')
        G_undirected = G_projected.to_undirected()
        G_edges_as_gdf = ox.graph_to_gdfs(G_undirected, nodes=False, edges=True)

        split_lines = split_linestring(G_edges_as_gdf)
        block_faces = list(polygonize(split_lines['geometry']))
        blocks = gpd.GeoDataFrame(geometry=block_faces)
        blocks = blocks.set_crs('EPSG:4326', allow_override=True)

        PolyInCity = sjoin(blocks, df_city, how='inner')
        PolyInCity.drop(columns=["index_right"], inplace=True)

        merged_polygons = gpd.GeoSeries(cascaded_union(PolyInCity["geometry"].values))
        merged_polygons = merged_polygons.set_crs('EPSG:4326', allow_override=True)

        rest = df_city.difference(merged_polygons)
        rest = gpd.GeoDataFrame(rest)
        rest = rest.rename(columns={0: 'geometry'}).set_geometry('geometry')

        rest_poly = explode(rest)
        rest_poly.reset_index(inplace=True)
        rest_poly.drop(columns=["level_0"], inplace=True)
        rest_poly.rename(columns={"level_1": "osm_id"}, inplace=True)

        new_df = PolyInCity.append(rest_poly, ignore_index=True)
        df_hc = new_df.to_crs('EPSG:5243')
        df_hc["centroid"] = df_hc.centroid
        coord = np.column_stack([df_hc["centroid"].x, df_hc["centroid"].y])

        # Think about automated resolution as threshold-> maybe depending of area size of boundary polygon
        threshold = get_LGU_threshold(df_city, 9)
        print(f"Threshold for nb of LGUs is {threshold}")
        th = get_hdscan_parameter(coord, threshold)

        model = AgglomerativeClustering(n_clusters=None, distance_threshold=th, affinity='euclidean')
        model.fit(coord)

        new_df['Cluster'] = model.labels_
        df_hc['Cluster'] = model.labels_

        merged_polys = []
        for idx in new_df["Cluster"].unique():
            tmp = new_df[new_df["Cluster"] == idx]
            polygons = tmp["geometry"].to_numpy()
            merged_polygon = gpd.GeoSeries(cascaded_union(polygons))
            merged_polys.append(merged_polygon[0])

        new_merged_polys = {'geometry': merged_polys}
        merged_polys_df = gpd.GeoDataFrame(new_merged_polys, crs="EPSG:4326")

        keep_df = merged_polys_df[merged_polys_df.geom_type == "Polygon"]
        To_explode = merged_polys_df[merged_polys_df.geom_type == "MultiPolygon"]
        explode_df = explode(To_explode)
        explode_df = explode_df.reset_index()
        explode_df.drop(columns=["level_0", "level_1"], inplace=True)
        overall = pd.concat([keep_df, explode_df])

        return overall
