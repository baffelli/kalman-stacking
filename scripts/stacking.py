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

from itertools import islice



def stack(inputs, outputs, threads, config, params, wildcards):

    ifgram = input.ifgrams[0]
    par = ifgram + "_par"
    current_if = gpf.gammaDataset(par, ifgram, dtype=gpf.type_mapping['FCOMPLEX'])
    avg_if = np.zeros(current_if.shape, dtype=np.complex64)
    # Compute reference phase
    ref_idx = cf.window_idx(current_if, [params.ridx, params.azidx], [params.ws, params.ws])
    dt_avg = 0
    for ifgram, master_slc, slave_slc in zip(input.ifgrams, zip(input.mli, input.mli[1::])):
        par = master_slc + '.par'
        current_if = gpf.gammaDataset(par, ifgram, dtype=gpf.type_mapping['FCOMPLEX'])
        dt = current_if.temporal_baseline[0]
        avg_if += current_if * current_if[ref_idx].conj()
        dt_avg += dt
    avg_if /= dt_avg
    rgb, norm, rest = vf.dismph(avg_if, peak=True)
    plt.imshow(rgb, aspect=1 / 10.0)
    plt.show()
    # TODO fix write_dataset so that it ssaves the correct type
    gpf.write_dataset(avg_if.astype(gpf.type_mapping['FCOMPLEX']), current_if.__dict__, output.avg_ifgram_par,
                      output.avg_ifgram)


stack(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)