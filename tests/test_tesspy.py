from unittest import TestCase
import tesspy.tesspy as tp
import geopandas as gpd
import pandas as pd
from shapely import wkt


class TestTessObj(TestCase):
    def test_quad_key(self):
        print("Test for string input")
        for res in [i for i in range(8, 17)]:
            qk = tp.TessObj.quadKey(self, "Frankfurt am Main", res)
            print(qk.shape)

        print("Test for GeoDataFrame")
        admin_polygon = tp.get_city_polygon("Frankfurt am Main")
        qk_GeoDF = tp.TessObj.quadKey(self, admin_polygon, 15)
        print(qk_GeoDF.shape)

        print("Test for Shapely Polygon")
        admin_polygon = tp.get_city_polygon("Frankfurt am Main")
        polygon = admin_polygon.geometry.values[0]
        qk_polygon = tp.TessObj.quadKey(self, polygon, 15)
        print(qk_polygon.shape)

    def test_hexagon(self):
        print("Test for string input")
        for res in [i for i in range(5, 10)]:
            h3 = tp.TessObj.hexagon(self, "Frankfurt am Main", res)
            print(h3.shape)

        print("Test for GeoDataFrame")
        admin_polygon = tp.get_city_polygon("Frankfurt am Main")
        h3_GeoDF = tp.TessObj.hexagon(self, admin_polygon, 8)
        print(h3_GeoDF.shape)

        print("Test for Shapely Polygon")
        admin_polygon = tp.get_city_polygon("Frankfurt am Main")
        polygon = admin_polygon.geometry.values[0]
        h3_polygon = tp.TessObj.hexagon(self, polygon, 8)
        print(h3_polygon.shape)

    # def test_voronoi(self):
    #     self.fail()

    # def test_cityblocks(self):
    #     self.fail()

    def test_admin_polygon(self):
        city = tp.get_city_polygon("Frankfurt am Main")
        print(city)

    def test_get_bbox(self):
        bbox = tp.get_bbox("Frankfurt am Main")
        print(bbox)

    def test_explode(self):
        qk = tp.TessObj.hexagon(self, "Frankfurt am Main", 14)
        exploded_qk = qk.explode()
        print(exploded_qk.head())

    def test_h3(self):
        area = tp.get_city_polygon("Frankfurt am Main")
        h3_tess = tp.get_h3(area, 9)
        print(h3_tess.shape)

    def test_adaptive_quadkey(self):
        admin_polygon = tp.get_city_polygon("Frankfurt am Main")
        test_data = pd.read_csv("FFM_clean_data.csv")
        test_data = gpd.GeoDataFrame(test_data, crs="EPSG:4326")
        test_data["geometry"] = test_data["geometry"].apply(wkt.loads)
        test_data.drop(columns=["Unnamed: 0"], inplace=True)
        test_data = test_data[["osmid", "geometry", "value"]]

        aqk = tp.TessObj.adaptive_squares(self, "Frankfurt am Main", test_data, 14)

        print(aqk.head())


    def test_cityblock(self):
        df = tp.TessObj.cityblocks(self, "Frankfurt am Main")
        print(df.shape)

    def test_voronoi(self):
        data = gpd.read_file("POI_data_FFM.geojson")
        df_vor = tp.TessObj.voronoi(self, "Frankfurt am Main", data, cluster_algo="hdbscan")
        print(df_vor.head(10))

    def test_voronoi_kmeans(self):
        data = gpd.read_file("POI_data_FFM.geojson")
        df_vor = tp.TessObj.voronoi(self, "Frankfurt am Main", data)
        print(df_vor.head(10))


    def test_square_polyfill(self):
        df_city = tp.get_city_polygon("Frankfurt am Main")
        gdf_square = tp.square_polyfill(df_city, 15)
        print(gdf_square.head())