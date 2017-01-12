feature_name = "Bisgletscher"
from osgeo import ogr
from osgeo import gdal
from pyrat.geo import geofun
from osgeo import gdalnumeric
import matplotlib.pyplot as plt
import gdal
import ogr
import osr
import cartopy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pyrat.geo.geofun as geo
import pyrat.visualization.visfun as vf

def overlay_displacement(inputs, outputs, threads, config, params, wildcards):
    #Load interferogram
    ifgram = gdal.Open(inputs.diff)
    #Load mask
    outline = ogr.Open(inputs.mask)
    outline_layer = outline.GetLayer()
    outline_layer.SetAttributeFilter("NAME = '{name}'".format(name=feature_name))
    #Load basemap
    bm = gdal.Open(inputs.basemap)
    #Get spatial reference from basemap
    ref = osr.SpatialReference()
    ref.ImportFromWkt(bm.GetProjection())
    #Create plot
    f, ax = plt.subplots(subplot_kw={'projection': ccrs.epsg(ref.GetAuthorityCode('PROJCS'))})
    #Display basemap
    ax.imshow(bm.ReadAsArray().transpose((1, 2, 0)), extent=geo.get_ds_extent(bm), origin='upper')
    rgb, map, a = vf.dismph(ifgram.ReadAsArray(), black_background=False, coherence=True, coherence_threshold=0.6)
    ax.imshow(rgb, extent=geo.get_ds_extent(ifgram), origin='upper')
    plt.show()


overlay_displacement(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)