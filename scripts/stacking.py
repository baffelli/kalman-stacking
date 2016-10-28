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

# import pyrat.stacks as stacks
# from itertools import islice


#Stacking interferograms with a moving window
def stack(inputs, outputs, threads, config, params, wildcards):
    #Class for coordinate transformation
    table = geofun.GeocodingTable(inputs.dem_par, inputs.lut)
    radar_coord = table.geo_coord_to_radar_coord(config['interferogram']['reference_coordinate'])
    print(radar_coord)
    #compute coordinate
    stack = []
    #Create the stack
    for ifgram_path in inputs.ifgrams:
        ifgram = diff.Interferogram(ifgram_path + '_par', ifgram_path, dtype=gpf.type_mapping['FCOMPLEX'])
        #reference
        ref_reg = cf.window_idx(ifgram, radar_coord, (16, 16))
        ifgram *= np.exp(-1j * np.mean(np.angle(ifgram[ref_reg])))
        stack.append(ifgram)
        # RGB, pal, norm = vf.dismph(ifgram)
        # plt.imshow(RGB)
        # plt.plot(radar_coord[1], radar_coord[0],'bo')
        # plt.show()

    ph_rate = np.zeros(ifgram.shape, dtype=np.complex64)
    tb = [x.temporal_baseline for x in stack]
    ph_rate = np.average(np.angle(np.dstack(stack)),weights=tb,axis=-1)
    av_cc = np.average(np.abs(np.dstack(stack)), weights=tb, axis=-1)
    total_t = np.sum(tb)
    # window_list = cf.moving_window(inputs.ifgrams, n=5)
    # stacked = []
    # for window_index, window_files in enumerate(window_list):
    #     current_stack = []
    #     for interferogram_index, interferogram_path in enumerate(window_files):
    #         ifgram = diff.Interferogram(interferogram_path + '_par', interferogram_path, dtype=gpf.type_mapping['FCOMPLEX'])
    #         print(ifgram.temporal_baseline)
    #         current_stack.append(ifgram)
    #     print("The current stack contains:")
    #     #Build datacube
    #     stack = np.dstack(current_stack)
    #     #Average interferograms
    #     avg_if = np.mean(stack,axis=-1)
    #     #Append to history
    #     stacked.append(avg_if)
    rgb, map, pal  = vf.dismph(av_cc * np.exp(1j*ph_rate))
    plt.imshow(rgb)
    plt.show()

stack(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)