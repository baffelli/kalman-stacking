import pyrat.fileutils.gpri_files as gpf
import numpy as np
from snakemake import shell
import json
import pyrat.geo.geofun as gf


def adf(input, output, threads, config, params, wildcards):
    import pyrat.fileutils.gpri_files as gpf
    wd = gpf.par_to_dict(input.mli_par)['range_samples']
    adf_cmd = "adf {input.ifgram} {output.int_sm} {output.cc} {wd} {params.adf_alpha} {params.adf_window}".format(
        wd=wd, output=output, input=input, params=params)
    shell(adf_cmd)


adf(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
    snakemake.wildcards)