import pyrat.diff.core as diff
import pyrat.diff.utils as utils
import numpy as np
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf

def cov(input, output, threads, config, params, wildcards):
    # #Create itab
    # itab=
    #Import the stack
    stack = diff.Stack(input.diff_pars, input.unw, input.mli_par)
    #Vector of temporal baselines
    # print(stack.dt)


    #Load one mli to display
    mli = gpf.gammaDataset(input.mli_par[-1], input.mli[-1])
    #Convert it to a matrix
    #slice stable_area
    stable_slice = (slice(0,1000),slice(0,None))
    pt = (854,38)
    pt_1 = (271, 71)
    stack_matrix = np.dstack(stack.stack)[stable_slice]
    #Reshape
    # stack_matrix = stack_matrix.reshape(stack_matrix.shape[0:2] + (stack_matrix.shape[-1]//5,5))
    # print(stack_matrix.shape)
    #Compute average ifgram
    avg = np.mean(stack_matrix, axis=-1)
    #Compute outer product
    stack_outer = cf.smooth(np.einsum('...i,...j->...ij', stack_matrix, stack_matrix.conj()),[5,5,1,1])
    outer_diag = 1/np.sqrt(np.diagonal(stack_outer,axis1=-2, axis2=-1))
    outer_coh= np.einsum('...j,...ij, ...i->...ij',outer_diag,stack_outer,outer_diag)
    names = ["{:%H %M}-{:%H %M}".format(s.master_time, s.slave_time) for s in stack.stack]
    #New figure
    gs = gridspec.GridSpec(4, 4)
    gs.update(hspace=0.5)
    f = plt.figure()
    ifgram_ax = f.add_subplot(gs[:,0:2])
    cov_ax =  f.add_subplot(gs[0:2,2::])
    cov_ax_1 = f.add_subplot(gs[2::, 2::])
    #Plot stacked data
    ifgram_ax.imshow(avg)
    ifgram_ax.plot(*(pt[::-1]), mfc='r', marker='o')
    ifgram_ax.plot(*(pt_1[::-1]), mfc='r', marker='o')
    #Plot covariance matrix
    cov_ax.matshow(np.abs(outer_coh[pt]))
    cov_ax_1.matshow(np.abs(outer_coh[pt_1]))
    cov_ax_1.set_xticklabels(names)
    cov_ax_1.set_yticklabels(names)


    plt.show()

cov(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
            snakemake.wildcards)


