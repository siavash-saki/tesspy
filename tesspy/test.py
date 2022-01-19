import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox

from tessellation import Tessellation
from poi_data import RoadData
# ffm_poly = ox.geocode_to_gdf('Frankfurt am Main')
# ffm = Tessellation(ffm_poly)
# print(ffm.area_gdf)

#ffm = Tessellation("Frankfurt am Main")
#print(ffm.area_gdf)
#road_network = RoadData(ffm.area_gdf)
#road_data = road_network.get_road_network()
#print(road_data.shape)


from time import sleep
from progressbar import progressbar


def gen_sleep(verbose=True):
    print('start')

    if verbose:
        print('verbose...')
        generator = progressbar(range(5))

    else:
        generator = range(10)

    print('before for loop..')
    for i in generator:
        sleep(1)

    print('done')


#gen_sleep()
