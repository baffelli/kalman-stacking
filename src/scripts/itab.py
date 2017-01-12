import pyrat.diff.intfun
import pyrat.ipt.core as ipt
import csv

def itab(inputs, outputs, threads, config, params, wildcards):
    nslc = len(inputs.slc_names)
    step = config['ptarg']['step']
    window = config['ptarg']['window']
    ref  =config['ptarg']['ref']
    stride = config['ptarg']['stride']
    with open(outputs.pitab, 'w+') as of:
        for line in pyrat.diff.intfun.itab(nslc, window, stride, step, ref):
            of.writelines(" ".join(map(str,line)) + " 1" + '\n')


itab(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)