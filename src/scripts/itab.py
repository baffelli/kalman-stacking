import pyrat.diff.intfun
import pyrat.ipt.core as ipt
import csv

def itab(inputs, outputs, threads, config, params, wildcards):
    nslc = len(inputs.slc_names)
    step = config['kalman']['step']
    window = config['kalman']['window']
    ref  =config['kalman']['ref']
    stride = config['kalman']['stride']
    with open(outputs.pitab, 'w+') as of:
        itab = pyrat.diff.intfun.itab(nslc, window, stride, step, ref)
        itab.tofile(of)


itab(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)