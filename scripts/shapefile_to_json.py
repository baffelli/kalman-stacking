from osgeo import ogr
import json
import pyrat.geo.geofun as gf

def shapefile_to_json(inputs, outputs, threads, config, params, wildcards):
    print('entered')
    #Create json feature collection
    feature_collection = {"type": "FeatureCollection",
                         "features": []
                         }
    #Create mapping table
    geocoding_table = gf.GeocodingTable(inputs.dem_par, inputs.lut)
    #Load shapefile
    shapefile = ogr.Open(inputs.reference_coord)
    points = shapefile.GetLayer(0)
    json_list = []
    for feature in points:
        geom = feature.GetGeometryRef()
        feature_json = json.loads(feature.ExportToJson())
        feature_json['id'] = feature.GetField('id')
        feature_json['properties']['radar_coordinates'] = list(map(int, geocoding_table.geo_coord_to_radar_coord([geom.GetX(),geom.GetY()])))
        feature_collection['features'].append(feature_json)
kls -ltr
with open(outputs.reference_coord, 'w+') as of:
        json.dump(feature_collection, of)
    #
    # x_min, x_max, y_min, y_max = outline_layer.GetExtent()
    #
    # #Set the out raster
    # goal_raster = gdal.GetDriverByName('GTiff').Create(outputs.mask, dem.RasterXSize, dem.RasterYSize, 1, gdal.GDT_Byte)
    # goal_raster.SetGeoTransform(gt)
    # goal_raster.SetProjection(outline_layer.GetSpatialRef().ExportToWkt())
    # band = goal_raster.GetRasterBand(1)
    # #write the shapefile
    # band.Fill(255)
    # gdal.RasterizeLayer(goal_raster,[1], outline_layer, burn_values=[0] )
    # array = band.ReadAsArray()

shapefile_to_json(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)