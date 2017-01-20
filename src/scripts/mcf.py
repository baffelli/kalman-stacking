import pyrat.fileutils.gpri_files as gpf
import numpy as np
from snakemake import shell
import json
import pyrat.geo.geofun as gf


def mcf(input, output, threads, config, params, wildcards):

    # # read location of reference point
    with open(input.reference_coord) as inputfile:
        ref_coord = gf.get_reference_coord(json.load(inputfile), 0)
        #        #get width of data
    wd = gpf.par_to_dict(input.mli_par)['range_samples']
    # unwrap data
    mcf_cmd = "mcf {input.ifgram} {input.cc} {input.cc_mask} {output.unw} {wd} {params.mode} - - - - 1 1 - {rinit} {azinit} 1 ".format(
        wd=wd, params=params, rinit=ref_coord[0], azinit=ref_coord[1], input=input, output=output)
    shell(mcf_cmd)



mcf(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
    snakemake.wildcards)