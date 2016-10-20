feature_name = "Bisgletscher"
from osgeo import ogr
from osgeo import gdal
from pyrat.geo import geofun
from osgeo import gdalnumeric
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
import numpy as np
import matplotlib.pyplot as plt
import pyrat.core.corefun as cf
import pyrat.diff.intfun as intfun

# import pyrat.stacks as stacks
# from itertools import islice



def stack(inputs, outputs, threads, config, params, wildcards):
    window_list = cf.moving_window(inputs.ifgrams, n=5)
    stacked = []
    for window_index, window_files in enumerate(window_list):
        current_stack = []
        for interferogram_index, interferogram_path in enumerate(window_files):
            ifgram = gpf.gammaDataset(interferogram_path + '.par', interferogram_path)
            current_stack.append(ifgram)
        #Build datacube
        stack = np.dstack(current_stack)
        #Average interferograms
        avg_if = np.mean(stack,axis=-1)
        #Append to history
        stacked.append(avg_if)
    rgb, map, pal  = vf.dismph(avg_if)
    plt.imshow(rgb)
    plt.show()

stack(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)