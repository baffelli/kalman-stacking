from pyrat.diff import kalman as ka
import pyrat.diff.intfun as intfun
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import scipy as sp
import numpy as np
import pyrat.fileutils.gpri_files as gpf
import os
import matplotlib.colors as cols

import pyrat.core.corefun as cf

import pyrat.geo.geofun as gf
import json

import pyrat.diff.core as ifgrams


import pyrat.visualization.visfun as vf

def kalman(input, output, threads, config, params, wildcards):
    #Load stack
    stacks = [ifgrams.Stack.fromfile(f) for f in input.k]
    ifgram_shape = z[-1].shape[0]
    x0 = np.tile([0.1,0.1], (ifgram_shape,1))

    #Initialize Kalman filter
    kf = ka.KalmanFilter(F=F[:,None,:,:],H=H[:,None,:,:])
    #Tune filter using EM algorithm
    kf.smooth(z)
    # kf.filter(z)






kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


