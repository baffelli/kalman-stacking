import numpy as np
import pyrat.fileutils.gpri_files as gpf
import snakemake
import os
# import pyrat.stacks as stacks
# from itertools import islice

# Stacking interferograms with a moving window
def interferogram(input, output, threads, config, params, wildcards):
    shell("create_offset {input.master_par} {input.slave_par} {output.int_par} - 2 1 0")
    shell(
        "SLC_intf {input.master} {input.slave} {input.master_par} {input.slave_par} {output.int_par} {output.ifgram} {params.rlks} {params.azlks} - - 0 0 1 - - - - -")
    ifgram_par = gpf.par_to_dict(output.int_par)
    # Load reference coord
    ref_coord = np.genfromtxt(input.reference_coord)
    # referencing
    ref_cmd = "cpx_max {{output.ifgram}} - {{output.ifgram_ref}} {wd} 0 {ridx} {azidx} {nr} {naz} - - - 1".format(
        ridx=ref_coord[0], azidx=ref_coord[1], nr=config['interferogram']['reference_region_size'][0],
        naz=config['interferogram']['reference_region_size'][1], wd=ifgram_par.interferogram_width)
    shell(ref_cmd)


interferogram(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
      snakemake.wildcards)
