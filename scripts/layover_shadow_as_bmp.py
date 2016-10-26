feature_name = "Bisgletscher"
from osgeo import ogr
from osgeo import gdal
from pyrat.geo import geofun
from osgeo import gdalnumeric
import pyrat.fileutils.gpri_files as gpf
import matplotlib.pyplot as plt
import numpy as  np
def as_bmp(inputs, outputs, threads, config, params, wildcards):
    ls, par = gpf.load_dataset(inputs.ref_mli_par, inputs.mask, dtype=gpf.type_mapping['UCHAR'])
    print(ls.mean())
    plt.imshow((ls == chr(16)).all())
    plt.show()



as_bmp(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)