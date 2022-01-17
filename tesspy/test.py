import numpy as np
import pandas as pd
import geopandas as gpd
import osmnx as ox

# from tessellation import Tessellation
# ffm_poly = ox.geocode_to_gdf('Frankfurt am Main')
# ffm = Tessellation(ffm_poly)
# print(ffm.area_gdf)

available=['amenity']

lst=['amenity', 'shop', 'building']

for i in lst:
    if i in available:
        print(i)

lst.remove('amenity')
print(lst)
