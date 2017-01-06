import pyrat.diff.utils
from pyrat.diff import kalman as ka
import pyrat.diff.intfun as intfun
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import scipy as sp
import numpy as np
import pyrat.fileutils.gpri_files as gpf
import os
import matplotlib.colors as cols

import pyrat.core.corefun as cf

import pyrat.geo.geofun as gf
import json

import pyrat.diff.core as ifgrams


import pyrat.visualization.visfun as vf

def create_itab(input, output, threads, config, params, wildcards):
    #Create itab
    itab = pyrat.diff.utils.Itab(config['kalman']['nstack'], window=config['ptarg']['window'], step=config['ptarg']['step'], stride=config['ptarg']['stride'], n_ref=config['ptarg']['ref'])
    itab.tofile(output.itab)


create_itab(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
            snakemake.wildcards)


