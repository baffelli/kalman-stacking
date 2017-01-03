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

def prepare_kalman(input, output, threads, config, params, wildcards):
    #load inputs
    dt = []
    master_times = []
    slave_times = []
    z = []
    #Create itab
    # itab = intfun.Itab
    #Read in the input stack
    stack = ifgrams.Stack(input.diff_pars, input.unw, input.mli_pars, input.itab)
    #get number of states
    nstates = config['kalman']['nstates']
    ifgram_shape = stack[0].shape
    #Reshape to npixels * size
    z = np.dstack(stack.stack).reshape((np.prod(ifgram_shape),) + (len(stack.stack),))
    stack_dt = stack[0].master_par.start_time - stack[-1].slave_par.start_time
    #Matrix for the stack transition
    F = intfun.F_model(stack_dt)
    #H matrix
    phase_factor = np.pi * 4 / gpf.lam(stack[0].master_par.radar_frequency)
    H_m = np.array([1,0]) * phase_factor
    H = stack.H_stack(intfun.F_model, H_m)
    #Save the matrices
    np.save(output.H, H)
    np.save(output.z, z)
    np.save(output.F, F)


prepare_kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


