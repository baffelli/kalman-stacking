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

import pyrat.geo.geofun as gf
import json

import pyrat.diff.core as ifgrams




def kalman(input, output, threads, config, params, wildcards):
    #Load reference position
    with open(input.reference_coord,'r') as infile:
        ref_pos = gf.get_reference_coord(json.load(infile))
    #load kalman state
    x = np.genfromtxt(input.x)
    P = np.genfromtxt(input.P)
    #load inputs
    dt = []
    master_times = []
    slave_times = []
    z = []
    #Create itab
    itab = np.genfromtxt(input.itab)
    #Read in the input
    for unw_par, unw in zip(input.diff_pars, input.unw):
        unw = ifgrams.Interferogram(unw_par, unw)
        z.append(unw[ref_pos])
        dt.append(unw.temporal_baseline)
        master_times.append(unw.master_mli_par.start_time)
        slave_times.append(unw.slave_mli_par.start_time)
    #Read mli parameters of first and last in stack
    stack_dt = master_times[0] - slave_times[-1]
    F = intfun.F_model(stack_dt)
    print(F)




kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


