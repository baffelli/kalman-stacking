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
    #Load kalman matrices
    z = np.array([np.load(cz) for cz in input.z])
    F = np.array([np.load(cf) for cf in input.F])
    H = np.array([np.load(ch) for ch in input.H])
    ifgram_shape = z[-1].shape[0]
    x0 = np.tile([0.1,0.1], (ifgram_shape,1))
    #Initialize Kalman filter
    R = np.tile(np.eye(z[-1].shape[-1]) * 1e-6, (z[-1].shape[0],1,1))
    Q = np.tile(np.eye(2) * 1e-3, (z[-1].shape[0],1,1))
    kf = ka.KalmanFilter(F=F,H=H, R=R, x0=x0, Q=Q)
    kf.filter(z)






kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


