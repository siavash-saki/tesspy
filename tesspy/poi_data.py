import numpy as np
import pandas as pd
import geopandas as gpd
import overpass
import warnings

def create_overpass_query_poi_data(area, poi_categories=['amenity']):
    """

    :param poi_categories: list
    :param area: GeoDataFrame
    :return: GeoDataFrame
    """
    # make the query area a bit larger
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        area = area.buffer(0.03).simplify(0.01)

    # make polygon string for OSM overpass query
    xy = np.array(area.iloc[0].exterior.coords)
    poly_str = ''
    for lat, lon in zip(xy[:, 1], xy[:, 0]):
        poly_str = poly_str + str(lat) + ' ' + str(lon) + ' '
    poly_str = poly_str.strip()

    # create query string for overpass
    query_string = ''
    for element in ['node', 'way']:
        for poi_category in poi_categories:
            query_string = query_string + f'{element}[{poi_category}](poly:"{poly_str}");'
    query_string = "[out:json];(" + query_string + ');out geom;'

    return query_string



def get_poi_data(query_string, timeout):

    api = overpass.API(timeout=timeout)
    resp = api._get_from_overpass(query).json()






    pass









def get_road_network(area):
    """
    :param area: GeoDataFrame
    :return: GeoDataFrame
    """
    pass