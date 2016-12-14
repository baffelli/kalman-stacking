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
    #load inputs
    dt = []
    master_times = []
    slave_times = []
    z = []
    #Create itab
    itab = intfun.Itab.fromfile(input.itab)
    #Read in the input stack
    stack = ifgrams.Stack(input.diff_pars, input.unw, input.itab)
    #load kalman state
    nstates = config['kalman']['nstates']
    ifgram_shape = stack[0].shape
    #Reshape to npixels * size
    z = np.array(stack.stack).reshape((np.prod(ifgram_shape),) + (len(stack.stack),))
    x = np.fromfile(input.x).reshape((np.prod(ifgram_shape),) + (nstates,))
    P = np.fromfile(input.P).reshape((np.prod(ifgram_shape),) + (nstates,nstates))
    #Read mli parameters of first and last in stack
    stack_dt = stack[0].master_par.start_time - stack[-1].slave_par.start_time
    #Master times
    master_times = [s.master_par.start_time for s in stack]
    #Matrix for the stack transition
    F = intfun.F_model(stack_dt)
    #H matrix
    phase_factor = np.pi * 4 / gpf.lam(stack[0].master_par.radar_frequency)
    H_m = np.array([1,0]) * phase_factor
    H = intfun.H_stack(intfun.F_model, H_m , itab,  master_times)
    #Covariance matrix
    R = np.eye(H.shape[0]) * 1e-3 * phase_factor
    Q = np.eye(2) * 1e-5
    print(z.shape)
    print(H.shape)
    filter = ka.KalmanFilter(2, H.shape[0], H=H[None,:,:], F=F[None,:,:], R=R[None,:,:], x0=x, Q=Q[None,:,:], P=P)
    #Prediction step
    filter.predict()
    #Update
    filter.update(z)
    filter.x.tofile(output.x)
    filter.P.tofile(output.P)




kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


