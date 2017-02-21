from osgeo import ogr
import json
import pyrat.geo.geofun as gf

def shapefile_to_json(inputs, outputs, threads, config, params, wildcards):
    #Create json feature collection
    feature_collection = {"type": "FeatureCollection",
                         "features": []
                         }
    #Create mapping table
    geocoding_table = gf.GeocodingTable(inputs.dem_par, inputs.lut, inputs.ref_mli, inputs.lut_inv)
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

    with open(outputs.reference_coord, 'w+') as of:
        json.dump(feature_collection, of)


shapefile_to_json(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)