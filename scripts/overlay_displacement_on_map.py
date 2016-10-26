feature_name = "Bisgletscher"
from osgeo import ogr
from osgeo import gdal
from pyrat.geo import geofun
from osgeo import gdalnumeric
import matplotlib.pyplot as plt

def overlay_displacement(inputs, outputs, threads, config, params, wildcards):
    #Load interferogram
    ifgram = gdal.Open(inputs.diff)
    #Load mask
    outline = ogr.Open(inputs.mask)
    outline_layer = outline.GetLayer()
    outline_layer.SetAttributeFilter("NAME = '{name}'".format(name=feature_name))


overlay_displacement(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)