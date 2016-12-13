from pyrat.diff import kalman as ka
import pyrat.diff.intfun as intfun
import numpy as np

def init_kalman(input, output, threads, config, params, wildcards):
    #Create itab and save it
    print(config)
    nstack = config['kalman']['nstack']
    itab = intfun.Itab(nstack, window=config['ptarg']['window'], step=config['ptarg']['step'], stride=config['ptarg']['stride'], n_ref=config['ptarg']['ref'])
    itab.tofile(output.itab)
    #Save initial filter state
    x = np.zeros(2)
    np.savetxt(output.x, x)
    #Save initial covariance
    P = np.eye(len(x))
    np.savetxt(output.P, P)



init_kalman(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)
