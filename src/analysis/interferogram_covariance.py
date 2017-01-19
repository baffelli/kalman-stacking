import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

import pyrat.core.corefun as cf
import pyrat.diff.core as diff
import pyrat.fileutils.gpri_files as gpf

import pyrat.geo.geofun as gf
import json

import scipy.misc as misc

def cov_to_image(cov_mat):
    image_shape = cov_mat.shape[0:2]
    cov_shape = cov_mat.shape[-2], cov_mat.shape[-1]
    new_shape = [image_shape[0] * cov_shape[0], image_shape[1] * cov_shape[1]]
    return cov_mat.reshape(new_shape)


def cov(input, output, threads, config, params, wildcards):
    plt.rcParams['font.size']=14
    # Import the stack
    stack = diff.Stack(input.diff_pars, input.unw, input.mli_par)
    # Load one mli to display
    mli = gpf.gammaDataset(input.mli_par[-1], input.mli[-1])
    #load reference
    with open(input.reference_coord) as inputfile:
       ref_coord =  gf.get_reference_coord(json.load(inputfile), 0)
    #load mask
    mask = misc.imread(input.mask).T
    # Convert it to a matrix
    # slice stable_area
    stable_slice = (slice(0, 1000), slice(0, None))
    # Slice MLI and interferograms
    mli = mli[stable_slice]
    stack_matrix = np.dstack(stack.stack)[stable_slice]
    mask = mask[stable_slice]
    # Compute average ifgram and phase variance
    avg = np.mean(stack_matrix, axis=-1)
    phase_var = np.var(stack_matrix, axis=-1)
    # detect the locations of maximum and minimum variance
    ref_idx = tuple(ref_coord)
    max_var_idx = (ref_coord[0]+50, ref_coord[1] + 40)
    # Compute outer product
    stack_outer = cf.smooth(np.einsum('...i,...j->...ij', stack_matrix, stack_matrix.conj()), [5, 5, 1, 1])
    outer_diag = 1 / np.sqrt(np.diagonal(stack_outer, axis1=-2, axis2=-1))
    outer_coh = np.einsum('...j,...ij, ...i->...ij', outer_diag, stack_outer, outer_diag)
    cov_im = cov_to_image(outer_coh[::5,::5])
    plt.imshow(np.abs(cov_im),vmin=0,vmax=1)
    plt.show()
    # l, w = np.linalg.eigh(stack_outer)
    # plt.imshow(np.abs(l[:,:,0]))
    # plt.show()
    # outer_l, outer_w = np.linalg.eigh(stack_outer)
    # Reshape
    outer_coh_im = cov_to_image(outer_coh)
    names = ["{:%H:%M}-{:%H:%M}".format(s.master_time, s.slave_time) for s in stack.stack]
    print(len(names))
    # New figure
    gs = gridspec.GridSpec(4, 5, width_ratios=[2, 2, 1, 1, 0.2])
    gs.update(hspace=0.9)
    f = plt.figure()
    # Plot mean and variance
    asp = 1 / 3
    var_ax = f.add_subplot(gs[::, 0:1])
    var_ax.imshow(phase_var, aspect=asp)
    mean_ax = f.add_subplot(gs[::, 1:2])
    mean_ax.imshow(avg, aspect=asp)
    mean_ax.imshow(mask, aspect=asp)
    mean_ax.plot(*(ref_idx[::-1]), marker='o', mfc='b')
    mean_ax.plot(*(max_var_idx[::-1]), marker='o', mfc='r')
    mean_ax.set_title('Average Phase')
    var_ax.set_title('Phase variance')
    # Add inset axes
    # Plot covariance matrices
    cov_ax_min = f.add_subplot(gs[0:2:, 2::4])
    cov_ax_min.imshow(np.abs(outer_coh[ref_idx]), vmin=0, vmax=1)
    cov_ax_min.set_title('Reference location')
    cov_ax_max = f.add_subplot(gs[2::, 2::4])
    cov_ax_max_mappable = cov_ax_max.imshow(np.abs(outer_coh[max_var_idx]), vmin=0, vmax=1)
    cov_ax_min.set_title('Test location')
    # Add colorbar
    cbar_ax = f.add_subplot(gs[::, -1])
    cbar = f.colorbar(cov_ax_max_mappable, cax=cbar_ax)
    cbar.set_label('Correlation coefficent')
    # Set ticks
    nstack =len(stack.stack)
    ticks_positions = range(nstack)
    for ax in [cov_ax_max, cov_ax_min]:
        ax.set_xticks(ticks_positions)
        ax.set_yticks(ticks_positions)
        ax.set_xticklabels(names, rotation='vertical')
        ax.set_yticklabels(names)
    # Connect covariance with point
    arrowprops = dict(facecolor='grey', arrowstyle='-')
    mean_ax.annotate('', xy=max_var_idx[::-1], xytext=(0, 0), xycoords=mean_ax.transData,
                     textcoords=cov_ax_max.transData, arrowprops=arrowprops)
    mean_ax.annotate('', xy=max_var_idx[::-1], xytext=(nstack, nstack), xycoords=mean_ax.transData,
                     textcoords=cov_ax_max.transData, arrowprops=arrowprops)

    mean_ax.annotate('', xy=ref_idx[::-1], xytext=(0, 0), xycoords=mean_ax.transData,
                     textcoords=cov_ax_min.transData, arrowprops=arrowprops)
    mean_ax.annotate('', xy=ref_idx[::-1], xytext=(nstack, nstack), xycoords=mean_ax.transData,
                     textcoords=cov_ax_min.transData, arrowprops=arrowprops)
    # ax.set_title('Interferogram Correlation')
    # cov_ax_1 = f.add_subplot(gs[2::, 2::])
    # #Plot stacked data
    # ifgram_ax.imshow(avg)
    # ifgram_ax.plot(*(pt[::-1]), mfc='r', marker='o')
    # ifgram_ax.plot(*(pt_1[::-1]), mfc='r', marker='o')
    # #Plot covariance matrix
    # cov_ax.matshow(np.abs(outer_coh[pt]))
    # cov_ax_1.matshow(np.abs(outer_coh[pt_1]))
    # cov_ax_1.set_xticklabels(names)
    # cov_ax_1.set_yticklabels(names)


    plt.show()
    f.savefig(output['cov_image'])


cov(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
    snakemake.wildcards)
