import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox


from tessellation import Tessellation
from tesspy.poi_data import RoadData
from tesspy.poi_data import POIdata
ffm = Tessellation("Kapstadt")
print(len(ffm.hexagons(8)))
#POI_ffm = POIdata(ffm_Polygon,  ["amenity", "building"], 60, False).get_poi_data()
#print(POI_ffm.shape)

#aqk = ffm.adaptive_squares(14, poi_categories=["amenity", "building"])
#print(aqk.shape)


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
