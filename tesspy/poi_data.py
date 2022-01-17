import numpy as np
import pandas as pd
from progressbar import progressbar
import overpass
import warnings


class POIdata:

    def __init__(self, area, poi_categories, timeout, verbose):
        self.area = area
        self.poi_categories = poi_categories
        self.timeout = timeout
        self.verbose = verbose

    def create_overpass_query_string(self):
        # make the query area a bit larger
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            area = self.area.buffer(0.03).simplify(0.01)

        # make polygon string for OSM overpass query
        xy = np.array(area.iloc[0].exterior.coords)
        poly_str = ''
        for lat, lon in zip(xy[:, 1], xy[:, 0]):
            poly_str = poly_str + str(lat) + ' ' + str(lon) + ' '
        poly_str = poly_str.strip()

        # create query string for overpass
        query_string = ''
        for element in ['node', 'way']:
            for poi_category in self.poi_categories:
                query_string = query_string + f'{element}[{poi_category}](poly:"{poly_str}");'
        query_string = "[out:json];(" + query_string + ');out geom;'

        return query_string

    def get_poi_data(self):

        query_string = self.create_overpass_query_string()

        if self.verbose:
            print('Getting data from OSM...')
        api = overpass.API(timeout=self.timeout)
        resp = api._get_from_overpass(query_string).json()

        lst_nodes = []
        lst_ways = []

        if self.verbose:
            print('Creating POI DataFrame...')
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


def get_road_network(area):
    """
    :param area: GeoDataFrame
    :return: GeoDataFrame
    """
    pass
