from sklearn.cluster import KMeans
import hdbscan
from scipy.spatial import Voronoi
from tessellation_functions import *
from poi_data import *
import numpy as np


# todo: return gdf with POI count

def get_city_polygon(city: str):
    """
    Gets the polygon of a city or an area

    Parameters
    ----------
    city : str
        city must be a name of a city, or an address of a region

    Returns
    --------
    df_city : GeoDataFrame
        GeoDataFrame containing the polygon of the city
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        df_city = ox.geocode_to_gdf(city)
    df_city = df_city[["osm_id", "geometry"]]
    df_city = df_city.rename({"osm_id": "osmid"})
    return df_city


def _check_input_geodataframe(gdf):
    """
    checks if the input gdf is in a correct format

    Parameters
    ----------
    gdf : GeoDataFrame

    Returns
    --------
    gdf : GeoDataFrame
    """
    if len(gdf) != 1:
        raise ValueError("GeoDataFrame must have only one geometry element")
    else:
        if not hasattr(gdf, "geometry"):
            raise TypeError("Geometry column missing in GeoDataFrame")
        else:
            if type(gdf['geometry'].iloc[0]) not in [Polygon, MultiPolygon]:
                raise TypeError("Geometry must be of type shapely polygon or multipolygon")
            else:
                if gdf.crs is None:
                    raise ValueError("GeoDataFrame must have a CRS")
                elif gdf.crs != 'epsg:4326':
                    return gdf.to_crs(epsg=4326)
                else:
                    return gdf


class Tessellation:
    """
    Creates a Tessellation object using a GeoDataFrame or a city name,
    which allows for different tessellation methods

    Parameters
    ----------
    area : GeoDataFrame or str
        GeoDataFrame must have a single shaply Polygon or MultiPolygon
        in geometry column and its CRS must be defined.
        str must be a name of a city, or an address of a region

    Examples
    --------
    >>> ffm= Tessellation('Frankfurt am Main')
    """

    def __init__(self, area):

        if type(area) == gpd.GeoDataFrame:
            self.area_gdf = _check_input_geodataframe(area)
        elif type(area) == str:
            self.area_gdf = get_city_polygon(area)
        else:
            raise TypeError("area must be in format: GeoDataFrame or string")

        self.poi_dataframe = pd.DataFrame()
        self.available_poi_categories = []
        self.road_network = pd.DataFrame()

    def _get_missing_poi_categories(self, poi_list):
        """
        Checks if the poi categories are already available in the object poi_dataframe
        Creates a list of missing categories, which should be downloaded from OSM

        Parameters
        ----------
        poi_list : list or 'all'
            A list of passed POI categories
            'all' means all the available POI categories

        Returns
        -------
        missing_poi_categories : list
            A list containing the POI categories which should be downloaded
        """
        if poi_list == 'all':
            poi_list = self.osm_primary_features()

        missing_poi_categories = []
        for poi_category in poi_list:
            if not hasattr(self.poi_dataframe, poi_category):
                missing_poi_categories.append(poi_category)
        return missing_poi_categories

    def squares(self, resolution: int):
        """
        Generate square grid laying over the area

        Parameters
        ----------
        resolution : int
            Specifies the size of squares
            A positive number between todo: complete

        Returns
        -------
        df_qk_squares : pandas.DataFrame
            Dataframe containing squares

        Examples
        --------
        todo: write examples
        """
        df_qk_squares = get_squares_polyfill(self.area_gdf, resolution)
        df_qk_squares = df_qk_squares.drop(columns=['osm_id', 'children_id'])

        return df_qk_squares

    def hexagons(self, resolution: int):
        """
        Generate hexagon grid laying over the area

        Parameters
        ----------
        resolution : int
            Specifies the size of hexagons
            A positive number between todo: complete

        Returns
        -------
        df_h3_hexagons : pandas.DataFrame
            Dataframe containing hexagons

        Examples
        --------
        todo: write examples
        """

        df_h3_hexagons = get_h3_hexagons(self.area_gdf, resolution)
        df_h3_hexagons = df_h3_hexagons.reset_index().rename(columns={'index': 'hex_id'})

        return df_h3_hexagons

    def adaptive_squares(self,
                         start_resolution: int,
                         poi_categories=["amenity", 'building'],
                         threshold=None,
                         timeout=60,
                         verbose=False):

        """
        Generate adaptive squares based on the input POI data.
        Squares are created at the start resolution. Each square is broken
        into four smaller squares while the number of its POI exceeds the
        threshold.
        POI categories should be a list of OSM primary map features.
        A complete list can be found on the OSM website:
        https://wiki.openstreetmap.org/wiki/Map_features

        Parameters
        ----------
        start_resolution : int
            Specifies the size of initial squares
            A positive number between todo: complete
        poi_categories : A list of OSM primary map features or 'all'
                         default=["amenity", 'building']
            'all' means all the available POI categories
            Possible values in the list: ['aerialway', 'aeroway',
            'amenity', 'barrier', 'boundary', 'building', 'craft',
            'emergency', 'geological', 'healthcare', 'highway',
            'historic', 'landuse', 'leisure', 'man_made', 'military',
            'natural', 'office', 'place', 'power', 'public_transport',
            'railway', 'route', 'shop', 'sport', 'telecom', 'tourism',
            'water', 'waterway']
        threshold : int, default=None
            Threshold for the number of POI in a single square. If square
            has more, it is divided into four squares. If None passed, the
            median number of POI per square in the initial level is used
            as threshold.
        timeout : int, default=60
            The TCP connection timeout for the request
        verbose : bool, default=False
            If True, print information while computing

        Returns
        -------
        df_adaptive_squares : pandas.DataFrame
            Dataframe containing adaptive squares

        Examples
        --------
        todo: write examples
        """

        # check if the categories are already available in poi_data
        # find the ones which are not available and should be downloaded from OSM
        missing_poi_categories = self._get_missing_poi_categories(poi_categories)

        # if there is any category that should be downloaded -->
        # create the query and run it
        if len(missing_poi_categories) > 0:
            poi_data_obj = POIdata(self.area_gdf, missing_poi_categories, timeout, verbose)
            poi_data_new = poi_data_obj.get_poi_data()
            # concat the data with the available poi_data dataframe of the Tessellation object
            self.poi_dataframe = pd.concat([self.poi_dataframe, poi_data_new]).fillna(False)
            self.poi_dataframe = self.poi_dataframe.reset_index(drop=True)

        # data, based on which, tessellation should be done
        tess_data = self.poi_dataframe[self.poi_dataframe[poi_categories].sum(axis=1) > 0]
        points_geom = tess_data[['center_longitude', 'center_latitude']] \
            .apply(lambda p: Point(p['center_longitude'], p['center_latitude']), axis=1)
        tess_data = gpd.GeoDataFrame(geometry=points_geom,
                                     data=tess_data[poi_categories],
                                     crs='EPSG:4326')
        poi_data_aqk = tess_data.rename(columns={"points_geom": "geometry"})

        df_aqk = get_squares_polyfill(self.area_gdf, start_resolution)
        aqk_count_df = count_poi(df_aqk, poi_data_aqk)

        if not threshold:
            threshold = int(np.median(aqk_count_df["count"].values))
            if verbose:
                print(f"Threshold={threshold}  ==> set as the median POI-count per square at the initial level")

        i = start_resolution
        while max(aqk_count_df["count"].values) > threshold:
            i += 1
            if verbose:
                print(f"Threshold exceeded! Squares are subdivided into resolution {i}")

            df_tmp = get_adaptive_squares(aqk_count_df, threshold)
            df_tmp.drop(columns=["count"], inplace=True)
            df_tmp2 = count_poi(df_tmp, poi_data_aqk)
            aqk_count_df = df_tmp2

        final_aqk = gpd.sjoin(aqk_count_df, self.area_gdf)
        final_aqk = final_aqk.drop(columns=['osm_id', 'children_id', 'index_right'])

        return final_aqk

    def voronoi(self,
                cluster_algo="k-means",
                poi_categories=["amenity", 'building'],
                timeout=60,
                n_polygons=100,
                min_cluster_size=15,
                verbose=False):
        """
        Generate Voronoi polygons based on the input POI data.
        POI categories should be a list of OSM primary map features.
        A complete list can be found on the OSM website:
        https://wiki.openstreetmap.org/wiki/Map_features

        Parameters
        ----------
        cluster_algo : {'k-means', 'hdbscan', None}, default='k-means'
            Algorithm for clustering the POI data before creating
            Voronoi generators. If None passed, POI data are
            directly used as generators.
        poi_categories : A list of OSM primary map features or 'all'
                         default=["amenity", 'building']
            'all' means all the available POI categories
            Possible values in the list: ['aerialway', 'aeroway',
            'amenity', 'barrier', 'boundary', 'building', 'craft',
            'emergency', 'geological', 'healthcare', 'highway',
            'historic', 'landuse', 'leisure', 'man_made', 'military',
            'natural', 'office', 'place', 'power', 'public_transport',
            'railway', 'route', 'shop', 'sport', 'telecom', 'tourism',
            'water', 'waterway']
        timeout : int, default=60
            The TCP connection timeout for the request
        n_polygons : int, default=100
            Only when cluster_algo="k-means", approximate number of
            polygons to be created. Positive number and less than
            the initial POI numbers
        min_cluster_size : int, default=15
            Only when cluster_algo="hdbscan", minimum cluster size
            Positive number
        verbose : bool, default=False
            If True, print information while computing

        Returns
        -------
        df_voronoi : pandas.DataFrame
            Dataframe containing Voronoi polygons

        Examples
        --------
        todo: write examples
        """

        # check if the categories are already available in poi_data
        # find the ones which are not available and should be downloaded from OSM
        missing_poi_categories = self._get_missing_poi_categories(poi_categories)

        if type(self.area_gdf) == MultiPolygon:
            queried_area = self.area_gdf.convex_hull
        else:
            queried_area = self.area_gdf

        # if there is any category that should be downloaded -->
        # create the query and run it
        if len(missing_poi_categories) > 0:
            poi_data_obj = POIdata(queried_area, missing_poi_categories, timeout, verbose)
            poi_data_new = poi_data_obj.get_poi_data()
            # concat the data with the available poi_data dataframe of the Tessellation object
            self.poi_dataframe = pd.concat([self.poi_dataframe, poi_data_new]).fillna(False)
            self.poi_dataframe = self.poi_dataframe.reset_index(drop=True)

        # data, based on which, tessellation should be done
        tess_data = self.poi_dataframe[self.poi_dataframe[poi_categories].sum(axis=1) > 0]
        data_locs = tess_data[['center_longitude', 'center_latitude']].values

        # create generators for Voronoi diagram
        if cluster_algo == "k-means":
            if verbose:
                print('K-Means Clustering...')
            clustering = KMeans(n_clusters=n_polygons).fit(data_locs)
            generators = [np.mean(data_locs[clustering.labels_ == label], axis=0) for label in range(n_polygons)]

        elif cluster_algo == "hdbscan":
            if verbose:
                print('HDBSCAN Clustering... This can take a while...')
            clustering = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, prediction_data=True).fit(data_locs)
            generators = [np.mean(data_locs[clustering.labels_ == label], axis=0) for label in
                          range(clustering.labels_.max() + 1)]

        elif cluster_algo is None:
            if len(tess_data["ClusterIDs"]) > 5000:
                raise ValueError("Too many generators for Voronoi diagram. Please select a clustering algorithm")
            else:
                generators = data_locs

        else:
            raise ValueError("Please use a clustering algorithm: k-means, hdbscan or None")

        # create Voronoi polygons
        if verbose:
            print('Creating Voronoi polygons...')
        voronoi_dia = Voronoi(generators)
        voronoi_poly = gpd.GeoDataFrame(geometry=[p for p in voronoi_polygons(voronoi_dia, 0.1)], crs="EPSG:4326")
        voronoi_poly = gpd.sjoin(voronoi_poly, self.area_gdf)
        vor_polygons = voronoi_poly.intersection(self.area_gdf.geometry.iloc[0])
        df_voronoi = gpd.GeoDataFrame(geometry=vor_polygons)

        return df_voronoi

    def city_blocks(self,
                    number_of_polygons=1000,
                    detail_deg=None,
                    split_roads=True,
                    verbose=False
                    ):
        """
        Create city bocks (tiles) using road data from the area. To collect road data and specify the highway types
        that should be included the custom filter can be used
        :param number_of_polygons: int, default = 1000
            targeted number of city blocks
        :param detail_deg: int, default = None
            define the number of the top (int) highway types to use in the custom filter
        :param split_roads: bool, default = True
            True if LineStrings should be split up such that each LineString contains exactly 2 Points
        :param verbose : bool, default=False
            If True, print information while computing
        :return: geopandas.GeoDataFrame
            GeoDataFrame with cityblock tiles
        """
        if type(self.area_gdf) == MultiPolygon:
            queried_area = self.area_gdf.convex_hull
        else:
            queried_area = self.area_gdf

        road_data = RoadData(queried_area, detail_deg, split_roads, verbose).get_road_network()
        if verbose:
            print(f"Road data is collected. Overall {len(road_data)} streets are included.")
            print("Creating initial city blocks using the road network data")

        blocks = create_blocks(road_data)

        polygons_in_area = gpd.sjoin(blocks, queried_area, how='inner')
        polygons_in_area.drop(columns=["index_right"], inplace=True)
        if verbose:
            print(f"Filtered out {len(blocks) - len(polygons_in_area)} polygons, that where not in the area.")

        print(f"{len(polygons_in_area)} polygons in area")
        rest_polygons = get_rest_polygon(polygons_in_area, queried_area)

        city_blocks = polygons_in_area.append(rest_polygons)
        city_blocks_hc = city_blocks.to_crs("EPSG:5243")
        city_blocks_hc["centroid"] = city_blocks_hc.centroid

        coordinates = np.column_stack([city_blocks_hc["centroid"].x, city_blocks_hc["centroid"].y])
        if verbose:
            print("Threshold for hierarchical clustering is computed.")
        th = get_hierarchical_clustering_parameter(coordinates, number_of_polygons)

        if th:
            if verbose:
                print(f"Distance threshold for clustering is {th}.")

            model = AgglomerativeClustering(n_clusters=None, distance_threshold=th, affinity='euclidean')
            model.fit(coordinates)
        else:
            if verbose:
                print(f"No threshold could match the inputs, so the nb of LGUs ({number_of_polygons}) is used.")

            model = AgglomerativeClustering(n_clusters=number_of_polygons, affinity='euclidean')
            model.fit(coordinates)

        city_blocks["Cluster"] = model.labels_
        city_blocks_hc["Cluster"] = model.labels_

        if verbose:
            print("Merging small city blocks...")
        merged_polys = []
        for idx in city_blocks["Cluster"].unique():
            tmp = city_blocks[city_blocks["Cluster"] == idx]
            polygons = tmp["geometry"].to_numpy()
            merged_polygon = gpd.GeoSeries(unary_union(polygons))
            merged_polys.append(merged_polygon[0])

        new_merged_polys = {'geometry': merged_polys}
        merged_polys_df = gpd.GeoDataFrame(new_merged_polys, crs="EPSG:4326")

        keep_df = merged_polys_df[merged_polys_df.geom_type == "Polygon"]
        To_explode = merged_polys_df[merged_polys_df.geom_type == "MultiPolygon"]
        explode_df = explode(To_explode)
        explode_df = explode_df.reset_index()
        explode_df.drop(columns=["level_0", "level_1"], inplace=True)
        final_city_blocks = pd.concat([keep_df, explode_df])

        return final_city_blocks

    def get_polygon(self):
        """
        returns the area polygon in GeoDataFrame format
        """
        return self.area_gdf

    def get_poi_data(self):
        """
        returns the POI data in GeoDataFrame format
        """
        return self.poi_dataframe

    def get_road_network(self):
        """
        returns the road network data in GeoDataFrame format
        """
        return self.road_network

    @staticmethod
    def osm_primary_features():
        """
        list of primary OSM features
        available at https://wiki.openstreetmap.org/wiki/Map_features

        Returns
        --------
        list
        """
        return POIdata.osm_primary_features()

    @staticmethod
    def osm_highway_types():
        """
        list of all highway types
        :return: list of highway types
        """
        return RoadData.osm_highway_types()
