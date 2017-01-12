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
import matplotlib.gridspec as gridspec
import datetime as dt


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
        ifgram = diff.Interferogram(ifgram_path + '_par', ifgram_path, dtype=gpf.type_mapping['FLOAT'])
        #reference
        ref_reg = cf.window_idx(ifgram, stable_cood, (16, 16))
        ifgram *= np.mean(ifgram[ref_reg])#reference to point
        stack.append(ifgram)
    ph_rate = np.zeros(ifgram.shape, dtype=np.float)
    # tb = [x.master_par.start_time - x.slave_par.start_time for x in stack]
    ws = 10
    moving_windows = cf.moving_window(stack, n=ws)
    mean = []
    variance = []
    time = []
    for window_index, window_files in enumerate(moving_windows):
        current_stack = []
        current_baseline = []
        current_time = []
        for interferogram_index, interferogram in enumerate(window_files):
            baseline = interferogram.slave_par.start_time - interferogram.master_par.start_time
            print("The current temporal baseline is {bl}".format(bl=baseline))
            current_stack.append(interferogram)
            current_baseline.append(baseline)
            master_time = dt.datetime.combine(interferogram.master_par.date, dt.time(0,0,0)) +  dt.timedelta(seconds=interferogram.master_par.start_time)
            current_time.append(master_time)
        current_stack = np.dstack(current_stack)
        #Compute mean and variance inside the window
        weights = np.array(current_baseline)
        weights[weights > 60*5] = 0
        current_mean = np.average(current_stack,axis=-1, weights=weights)
        mean.append(current_mean)
        variance.append(np.average((current_stack - current_mean[:,:,None])**2, weights=weights,axis=-1))
        time.append(current_time[len(current_time)//2])
        # print(current_stack.shape)
    variance = np.dstack(variance)
    mean = np.dstack(mean)
    point_variance = variance[stable_cood]
    point_mean = mean[stable_cood]
    gs = gridspec.GridSpec(2, 2)
    plot_ax = plt.subplot(gs[0,::])
    plot_ax.plot(time, point_mean)
    plt.fill_between(time, point_mean - point_variance**0.5, point_mean + point_variance**0.5, alpha=0.5)
    disp_ax = plt.subplot(gs[1::,0])
    # rgb, norm, pal = vf.dismph(mean[:,:,0])
    disp_ax.imshow(variance[:,:,5])
    disp_ax.plot(stable_cood[1], stable_cood[0], 'ro')
    plt.show()

stack(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)