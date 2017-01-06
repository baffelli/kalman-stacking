import pyrat.diff.utils
from pyrat.diff import kalman as ka
import pyrat.diff.intfun as intfun
import pyrat.fileutils.gpri_files as gpf
import numpy as np

def init_kalman(input, output, threads, config, params, wildcards):
    #Create itab and save it
    nstack = config['kalman']['nstack']
    nstates = config['kalman']['nstates']
    itab = pyrat.diff.utils.Itab(nstack, window=config['ptarg']['window'], step=config['ptarg']['step'], stride=config['ptarg']['stride'], n_ref=config['ptarg']['ref'])
    itab.tofile(output.itab)
    #Load mli par
    mli_par = gpf.par_to_dict(input.mli_par)
    ifgram_shape = (mli_par.range_samples, mli_par.azimuth_lines)
    #Save initial filter state
    x = np.zeros((np.prod(ifgram_shape),) + (nstates,))
    x[:,1] = 1e-6
    x.tofile(output.x)
    #Save initial covariance
    P = np.tile(np.eye(nstates) * 1e-3, (np.prod(ifgram_shape),) + (1,1))
    P.tofile(output.P)



init_kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)
