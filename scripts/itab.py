
import pyrat.ipt.core as ipt
import csv

def itab(inputs, outputs, threads, config, params, wildcards):
    if config['ptarg']['method'] == 'single_reference':
        sl = config['ptarg']['reference_number']
        master_name = inputs.slc_names[sl]
        idx_master = sl
        master_it = ([sl], [master_name])
    else:


    with open(outputs.pitab, 'w+') as of:
        counter = 0
        for idx_master, master in master_it:
            for idx_slave, slave in enumerate(inputs.slc_names):
                if idx_slave > idx_master:
                    of.write("{idx_master} {idx_slave} {counter} 1  {master} {slave}\n".format(idx_master=idx_master, idx_slave=idx_slave,
                                                                                         counter=counter, master=master, slave=slave)


itab(snakemake.input, snakemake.output, snakemake.threads, snakemake.config, snakemake.params,
              snakemake.wildcards)