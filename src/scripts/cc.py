import pyrat.fileutils.gpri_files as gpf
import numpy as np
from snakemake import shell
import json
import pyrat.geo.geofun as gf


def cc(input, output, threads, config, params, wildcards):
    wd = gpf.par_to_dict(input.mli1_par)['range_samples']
    cc_cmd = "cc_wave {input.ifgram} {input.mli1} {input.mli2} {output.cc} {wd} {params.rlks} {params.azlks} 0".format(
        wd=wd, input=input, output=output, params=params)
    shell(cc_cmd)


cc(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
    snakemake.wildcards)