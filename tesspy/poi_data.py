import numpy as np
import pandas as pd
import overpass
import warnings
import osmnx as ox
from tessellation_functions import split_linestring


class POIdata:
    """
    This class creates a query for the investigated area and POI categories.
    The query is sent to osm using overpass API and the data is retrieved.

    Parameters
    ----------
    area : GeoDataFrame or str
        GeoDataFrame must have a single shapely Polygon or MultiPolygon
        in geometry column and its CRS must be defined.
        str must be a name of a city, or an address of a region

    poi_categories : A list of OSM primary map features or 'all'
    timeout : int
        The TCP connection timeout for the overpass request
    verbose : bool
        If True, print information while computing
    """

    def __init__(self, area, poi_categories, timeout, verbose):
        self.area = area
        self.poi_categories = poi_categories
        self.timeout = timeout
        self.verbose = verbose

    @staticmethod
    def osm_primary_features():
        """
        list of primary OSM features
        available at https://wiki.openstreetmap.org/wiki/Map_features

        Returns
        --------
        osm_primary_features_lst : list
        """
        osm_primary_features_lst = ['aerialway',
                                    'aeroway',
                                    'amenity',
                                    'barrier',
                                    'boundary',
                                    'building',
                                    'craft',
                                    'emergency',
                                    'geological',
                                    'healthcare',
                                    'highway',
                                    'historic',
                                    'landuse',
                                    'leisure',
                                    'man_made',
                                    'military',
                                    'natural',
                                    'office',
                                    'place',
                                    'power',
                                    'public_transport',
                                    'railway',
                                    'route',
                                    'shop',
                                    'sport',
                                    'telecom',
                                    'tourism',
                                    'water',
                                    'waterway']

        return osm_primary_features_lst


    def create_overpass_query_string(self):
        """
        creates the query string to be passed to overpass

        Returns
        --------
        query_string : str
        """
        # make the query area a bit larger
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            area = self.area.buffer(0.008).simplify(0.005)

        # make polygon string for OSM overpass query
        xy = np.array(area.iloc[0].exterior.coords)
        poly_str = ''
        for lat, lon in zip(xy[:, 1], xy[:, 0]):
            poly_str = poly_str + str(lat) + ' ' + str(lon) + ' '
        poly_str = poly_str.strip()

        # if poi not in primary --> error
        # todo: is this necessary?
        for poi_category in self.poi_categories:
            if poi_category not in self.osm_primary_features():
                raise ValueError(f'{poi_category} is not a valid POI primary category. See a list of OSM primary '
                                 f'features with Tessellation.osm_primary_features()')

        # create query string for overpass
        query_string = ''
        for element in ['node', 'way']:
            for poi_category in self.poi_categories:
                query_string = query_string + f'{element}[{poi_category}](poly:"{poly_str}");'
        query_string = "[out:json];(" + query_string + ');out geom;'

        return query_string

    def get_poi_data(self):
        """
        sends the query to osm using the overpass API and gets the data

        Returns
        --------
        poi_df : pandas.DataFrame
            A dataframe containing the POI, POI type, and coordinates
        """

        query_string = self.create_overpass_query_string()

        if self.verbose:
            print('Getting data from OSM...')
        api = overpass.API(timeout=self.timeout)
        resp = api._get_from_overpass(query_string).json()

        if self.verbose:
            print('Creating POI DataFrame...')

        lst_nodes = []
        lst_ways = []

        # if self.verbose:
        #     generator = tqdm(resp['elements'])
        # else:
        #     generator = resp['elements']

        generator = resp['elements']

        for item in generator:

            for cat in self.poi_categories:
                if cat in item['tags'].keys():
                    item[cat] = True
            if item['type'] == 'node':
                lst_nodes.append(item)
            elif item['type'] == 'way':
                item['center_latitude'] = np.mean([point['lat'] for point in item['geometry']])
                item['center_longitude'] = np.mean([point['lon'] for point in item['geometry']])
                lst_ways.append(item)
            else:
                continue

        if self.verbose:
            print('Cleaning POI DataFrame...')

        nodes_df = pd.DataFrame(lst_nodes)
        ways_df = pd.DataFrame(lst_ways)

        nodes_df['geometry'] = nodes_df[['lon', 'lat']].apply(lambda p: [{'lat': p['lat'], 'lon': p['lon']}], axis=1)
        nodes_df = nodes_df.rename(columns={'lat': 'center_latitude', 'lon': 'center_longitude'})
        nodes_df = nodes_df.drop(columns=['id'])
        ways_df = ways_df.drop(columns=['id', 'bounds', 'nodes'])

        poi_df = pd.concat([ways_df, nodes_df]).fillna(False)

        first_cols = ['type', 'geometry', 'tags', 'center_latitude', 'center_longitude']
        second_cols = sorted(poi_df.columns.drop(first_cols))
        poi_df = poi_df[first_cols + second_cols]

        return poi_df


class RoadData:
    """
        This class creates a custom filter for the investigated area and highway types.
        The custom filter is used to collect road data using osmnx.

        Parameters
        ----------
        area : GeoDataFrame or str
            GeoDataFrame must have a single shapely Polygon or MultiPolygon
            in geometry column and its CRS must be defined.
            str must be a name of a city, or an address of a region

        detail_deg : int
            integer to select the top (int) highway types
        split_roads : bool
            decide if LineStrings will be split up such that each LineString contains exactly 2 Points
        verbose : bool
            If True, print information while computing
    """

    def __init__(self, area, detail_deg=None, split_roads=True, verbose=False):
        self.detail_deg = detail_deg
        self.area = area
        self.verbose = verbose
        self.split_roads = split_roads

    @staticmethod
    def osm_highway_types():
        """
        list of OSM highway types
        available at https://wiki.openstreetmap.org/wiki/Key:highway

        Returns
        --------
        osm_highways_lst : list
        """
        osm_highways_lst = ['motorway',
                            'trunk',
                            'primary',
                            'secondary',
                            'tertiary',
                            'residential',
                            'unclassified',
                            'motorway_link',
                            'trunk_link',
                            'primary_link',
                            'secondary_link',
                            'living_street',
                            'pedestrian',
                            'track',
                            'bus_guideway',
                            'footway',
                            'path',
                            'service',
                            'cycleway']
        return osm_highways_lst

    def create_custom_filter(self):
        """
        Creates custom filter that is used to collect road data with osmnx. The detail_deg is used to use only the top
        highway types.
        :return: string with highway types that is used to collect data using osmnx
        """
        if self.detail_deg is None:
            if self.verbose:
                print("Creating custom filter for all highway types")
            highwaytypes = self.osm_highway_types()
        elif type(self.detail_deg) is int:
            if self.verbose:
                print("Creating custom filter for the specified detail degree")
            highwaytypes = self.osm_highway_types()[:self.detail_deg]
        else:
            raise ValueError("Please insert a valid detail degree: None or int")

        query = ''
        for types in highwaytypes[:-1]:
            query += f'{types}|'
        query += highwaytypes[-1]
        custom_filter = f"['highway'~'{query}']"

        if self.verbose:
            print(f"Created custom filter is {custom_filter}")
        return custom_filter

    def get_road_network(self):
        """
        Collects the road network data based on the defined custom filter. The inital roaddata is based on graphs
        and transformed into a GeoDataFrame. Within the RoadData Class the user can devide to split the data such
        that each LineString (representing a streetsegment) contains only two points and not multiple points.
        :return: GeoDataFrame containing road network
        """
        cf = self.create_custom_filter()

        if self.verbose:
            print("Collection street network data")
        G = ox.graph_from_polygon(self.area.boundary.convex_hull.values[0], custom_filter=cf)
        G_projected = ox.project_graph(G, to_crs='epsg:4326')
        G_undirected = G_projected.to_undirected()
        G_edges_as_gdf = ox.graph_to_gdfs(G_undirected, nodes=False, edges=True)
        if self.verbose:
            print("Splitting the linestring, such that each linestring has exactly 2 points.")

        if self.split_roads:
            road_network = split_linestring(G_edges_as_gdf)
        else:
            road_network = G_edges_as_gdf
        if self.verbose:
            print(f"Collected data has {len(road_network)} street segments")

        return road_network
