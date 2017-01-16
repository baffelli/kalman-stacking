import pyrat.diff.core as diff
import pyrat.diff.utils as utils
import numpy as np
import pyrat.core.corefun as cf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import pyrat.fileutils.gpri_files as gpf
import pyrat.visualization.visfun as vf

def create_stack(input, output, threads, config, params, wildcards):
    # Import the stack
    itab = utils.Itab(len(input.unw), **config['kalman'])
    stack = diff.Stack(input.diff_pars, input.unw, input.mli_par)
    stack.itab = itab
    stack.tofile(output.stack)

create_stack(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
            snakemake.wildcards)