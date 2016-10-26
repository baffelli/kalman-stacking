feature_name = "Bisgletscher"
from osgeo import ogr
from osgeo import gdal
from pyrat.geo import geofun
from osgeo import gdalnumeric
import pyrat.fileutils.gpri_files as gpf
import matplotlib.pyplot as plt

def as_bmp(inputs, outputs, threads, config, params, wildcards):
    mask, mask_data = gpf.load_dataset(inputs.ref_mli_par, inputs.mask, dtype=gpf.type_mapping['UCHAR'])



as_bmp(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)