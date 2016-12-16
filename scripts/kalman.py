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


import pyrat.visualization.visfun as vf

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
    stack = ifgrams.Stack(input.diff_pars, input.unw, input.mli_pars, input.itab)
    #load kalman state
    nstates = config['kalman']['nstates']
    ifgram_shape = stack[0].shape
    #Reshape to npixels * size
    z = np.dstack(stack.stack).reshape((np.prod(ifgram_shape),) + (len(stack.stack),))
    #Sample covariance of the atmosphere
    R_samp = np.einsum('...i,...j->...ij',z,z.conj()).reshape(ifgram_shape +  (len(stack.stack),len(stack.stack)))

    plt.imshow(R_samp[:,:,2,0])

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
    H = stack.H_stack(intfun.F_model, H_m)
    #Covariance matrix
    R = np.eye(H.shape[0]) * 2 * phase_factor + np.ones((H.shape[0],H.shape[0])) * 0.2
    Q = np.eye(2) * 1e-6
    filter = ka.KalmanFilter(2, H.shape[0], H=H[None,:,:], F=F[None,:,:], R=R[None,:,:], x0=x, Q=Q[None,:,:], P=P)
    #Prediction step
    filter.predict()
    #Update
    filter.update(z)
    filter.x.tofile(output.x)
    filter.P.tofile(output.P)
    #Display
    seconds_to_day = 24 * 60 * 60
    f, (ax_d, ax_v) = plt.subplots(1,2)
    asp = 1/3
    # ax_z.imshow(z.reshape(ifgram_shape +  (len(stack.stack),))[:,:,0], aspect=asp)
    # ax_z.set_title('Unwrapped interferogram')
    displacement = x.reshape(ifgram_shape + (nstates,))[:,:,0]
    displacement_variance = P.reshape(ifgram_shape + (nstates,nstates))[:,:,0,0]
    print(displacement_variance)
    rgb = vf.disp_value_and_variance(displacement, displacement_variance, var_tresh=0.2)
    cax = ax_d.imshow(displacement, vmin=-2e-1, vmax=2e-1, aspect=asp,  cmap='RdBu_r')
    f.colorbar(cax, ax=ax_d)
    ax_d.set_title('Cumulative displacement [m]')
    velocity = x.reshape(ifgram_shape + (nstates,))[:,:,1] * seconds_to_day
    velocity_variance = P.reshape(ifgram_shape + (nstates,nstates))[:,:,1,1] * seconds_to_day**2
    rgb = vf.disp_value_and_variance(velocity, velocity_variance, var_tresh=0.2, vmin=-0.5, vmax=0.5)
    cax = ax_v.imshow(velocity, aspect=asp, vmin=-2, vmax=2, cmap='RdBu_r')
    f.colorbar(cax, ax=ax_v)
    ax_v.set_title('Displacement velocity [m/day]')
    #Times of first and last slc
    first_time =  stack[0].master_time
    last_time = stack[-1].slave_time
    title = "Stack on {} between {} and {}".format(first_time.date(), first_time.time(), last_time.time())
    f.suptitle(title)
    # plt.show()
    f.savefig(output.fig, dpi=200)




kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)


