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
    print(P)
    #load inputs
    dt = []
    master_times = []
    slave_times = []
    z = []
    #Create itab
    itab = intfun.Itab.fromfile(input.itab)
    #Read in the input
    for unw_par, unw in zip(input.diff_pars, input.unw):
        unw = ifgrams.Interferogram(unw_par, unw)
        z.append(unw[ref_pos[0] + 30, ref_pos[1] + 30])
        dt.append(unw.temporal_baseline)
        master_times.append(unw.master_par.start_time)
        slave_times.append(unw.slave_par.start_time)
    print(z)
    #Read mli parameters of first and last in stack
    stack_dt = master_times[0] - slave_times[-1]
    #Matrix for the stack transition
    F = intfun.F_model(stack_dt)
    #H matrix
    phase_factor = np.pi * 4 / gpf.lam(unw.master_par.radar_frequency)
    H_m = np.array([1,0]) * phase_factor
    H = intfun.H_stack(intfun.F_model, H_m , itab,  master_times)
    #Covariance matrix
    R = np.eye(H.shape[0]) * 1e-3 * phase_factor
    Q = np.eye(2) * 1e-5
    filter = ka.KalmanFilter(2, H.shape[0], H=H, F=F, R=R, x0=x, Q=Q, P=P)
    #Prediction step
    filter.predict()
    #Update
    filter.update(z)
    x1 = np.savetxt(output.x, filter.x)
    P1 = np.savetxt( output.P, filter.P)




kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


