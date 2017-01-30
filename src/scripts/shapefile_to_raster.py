feature_name = "Bisgletscher"
from osgeo import ogr
from osgeo import gdal

def shapefile_to_raster(inputs, outputs, threads, config, params, wildcards):

    #open DEM
    dem = gdal.Open(inputs.dem)
    gt = dem.GetGeoTransform()
    pixel_width = gt[1]
    pixepixel_height = gt[5]

    #Select only desired feature
    shapefile_path = '/vsizip/' + inputs.mask
    outline = ogr.Open(shapefile_path)
    print(outline)
    outline_layer = outline.GetLayer()
    outline_layer.SetAttributeFilter("NAME = '{name}'".format(name=params.feature_name))

    x_min, x_max, y_min, y_max = outline_layer.GetExtent()

    #Set the out raster
    goal_raster = gdal.GetDriverByName('GTiff').Create(outputs.mask, dem.RasterXSize, dem.RasterYSize, 1, gdal.GDT_Byte)
    goal_raster.SetGeoTransform(gt)
    goal_raster.SetProjection(outline_layer.GetSpatialRef().ExportToWkt())
    band = goal_raster.GetRasterBand(1)
    #write the shapefile
    band.Fill(255)
    gdal.RasterizeLayer(goal_raster,[1], outline_layer, burn_values=[0] )
    array = band.ReadAsArray()

shapefile_to_raster(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)