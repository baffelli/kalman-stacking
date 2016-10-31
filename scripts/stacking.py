from osgeo import ogr
from osgeo import gdal
from pyrat.geo import geofun
from osgeo import gdalnumeric
import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf
import numpy as np
import matplotlib.pyplot as plt
import pyrat.core.corefun as cf
import pyrat.diff.core as diff
import pyrat.diff.intfun as intfun
import pandas as pd


# import pyrat.stacks as stacks
# from itertools import islice

import scipy as sp
from scipy import signal as sig

#Stacking interferograms with a moving window
def stack(inputs, outputs, threads, config, params, wildcards):
    #Class for coordinate transformation
    table = geofun.GeocodingTable(inputs.dem_par, inputs.lut)
    stable_cood = table.geo_coord_to_radar_coord(config['interferogram']['reference_coordinate'])
    glacier_coord = table.geo_coord_to_radar_coord([	623006.4, 106231.6])
    #compute coordinate
    stack = []
    #Create the stack
    for ifgram_path in inputs.ifgrams:
        ifgram = diff.Interferogram(ifgram_path + '_par', ifgram_path, dtype=gpf.type_mapping['FCOMPLEX'])
        #reference
        ref_reg = cf.window_idx(ifgram, stable_cood, (16, 16))
        ifgram *= np.exp(-1j * np.mean(np.angle(ifgram[ref_reg])))#reference to point
        stack.append(ifgram)


    ph_rate = np.zeros(ifgram.shape, dtype=np.complex64)
    tb = [x.master_par.start_time - x.slave_par.start_time for x in stack]
    moving_windows = cf.moving_window(stack, n=15)
    mean = []
    variance = []
    time = []
    for window_index, window_files in enumerate(moving_windows):
        current_stack = []
        current_baseline = []
        current_time = []
        for interferogram_index, interferogram in enumerate(window_files):
            baseline = interferogram.master_par.start_time - interferogram.slave_par.start_time
            current_stack.append(interferogram)
            current_baseline.append(baseline)
            current_time.append( interferogram.master_par.start_time)
        current_stack = np.dstack(current_stack)
        #Compute mean and variance inside the window
        mean.append(np.average(current_stack,axis=-1))
        variance.append(np.var(np.angle(current_stack),axis=-1))
        time.append(np.mean(current_time))
        # print(current_stack.shape)
    variance = np.dstack(variance)
    mean = np.dstack(mean)
    print(mean.shape)
    plt.plot(time, np.angle(mean[stable_cood]))
    plt.show()

stack(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)