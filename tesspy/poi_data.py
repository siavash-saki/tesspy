import numpy as np
import pandas as pd
from progressbar import progressbar
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

    @staticmethod
    def osm_highway_types():
        """
        list of primary OSM features
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

    def create_custom_filter(self, detail_deg=None):

        if detail_deg is None:
            highwaytypes = self.osm_highway_types()
        elif type(detail_deg) is int:
            highwaytypes = self.osm_highway_types()[:detail_deg]
        else:
            raise ValueError("Please insert a valid detail degree: None or int")

        query = ''
        for types in highwaytypes[:-1]:
            query += f'{types}|'
        query += highwaytypes[-1]
        custom_filter = f"['highway'~'{query}']"
        return custom_filter

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

        if self.verbose:
            generator = progressbar(resp['elements'])
        else:
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

    def get_road_network(self, detail_deg=None):
        cf = self.create_custom_filter(detail_deg)
        box = self.area.geometry.values[0].minimum_rotated_rectangle
        x, y = box.exterior.coords.xy
        bbox = [max(y), min(y), max(x), min(x)]

        G = ox.graph_from_bbox(bbox[0], bbox[1], bbox[2], bbox[3], custom_filter=cf)
        G_projected = ox.project_graph(G, to_crs='epsg:4326')
        G_undirected = G_projected.to_undirected()
        G_edges_as_gdf = ox.graph_to_gdfs(G_undirected, nodes=False, edges=True)

        road_network = split_linestring(G_edges_as_gdf)
        return road_network
