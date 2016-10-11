import pyrat
from snakemake.remote.SFTP import RemoteProvider, RemoteObject
import datetime as dt
import csv
provider = RemoteProvider(port=22,username="baffelli", private_key="/home/baffelli/.ssh/id_rsa")

#tar = provider.glob_wildcards("ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/raw/{date,20150803}_{time,06\d+}.tar")

configfile: './bisgletscher.json'
slc_regex = "(\d{8})_(\d{6})((_[A-B]{3}[l-u]))?"
dt_regex = "(\d{8})_(\d{6})"
dtfmt = "%Y%m%d_%H%M%S"


def select_date_range(string_dates, date_start, date_end):
    dates = [ dt.datetime.strptime(x, dtfmt) for x in string_dates]
    dt_start = dt.datetime.strptime(date_start, dtfmt)
    dt_end = dt.datetime.strptime(date_end, dtfmt)
    valid_dates = [date.strftime(dtfmt) for date in dates if dt_start < date < dt_end]
    print(valid_dates)
    return valid_dates

def valid_dates_on_server(start_dt, stop_dt):
    with open('list_of_slcs.csv', 'r') as f:
        reader = csv.reader(f)
        all_dates = list(reader)
    valid_dates = select_date_range(all_dates[0], start_dt, stop_dt)
    return valid_dates

def all_ifgrams(wildcards):
    import itertools as iter
    valid_dates = valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
    ifgrams = expand(expand('diff/{{vd1}}_{chan}_{{vd2}}_{chan}.int', chan=wildcards.chan),zip, vd1=valid_dates, vd2=valid_dates[1::])
    print(ifgrams)
    return ifgrams





include: pyrat.rules['raw_to_slc']


rule all:
    input:
       'outputs/20150803_060019_20150803_063519_AAAl.int',
       'list_of_slcs.csv'





#Do not need to perform RC if the data comes from the server
ruleorder: untar_and_copy > range_compression
rule untar_and_copy:
        input:
            tar = provider.remote("ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/raw/{datetime}.tar", static=True),
        output:
            slc = 'slc_chan/{datetime}_{chan}.slc',
            slc_par = 'slc_chan/{datetime}_{chan}.slc.par'
        run:
            shell("tar -xvf {input.tar} --directory ./slc_chan --strip=3 data/Dom/slc/{wildcards.datetime}_{wildcards.chan}.slc")
            shell("tar -xvf {input.tar} --directory ./slc_chan --strip=3 data/Dom/slc/{wildcards.datetime}_{wildcards.chan}.slc.par ")

#Get a list of all slcs on the server
rule get_list_of_all_slcs:
    output:
        'list_of_slcs.csv'
    run:
        slcs = provider.glob_wildcards("ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/raw/{slcname}.tar")
        with open(output[0], 'w+') as of:
            writer = csv.writer(of)
            writer.writerows(slcs)


rule ifgram:
    input:
        master = 'slc_desq/{mastername}.slc_dec',
        master_par = 'slc_desq/{mastername}.slc_dec.par',
        slave = 'slc_desq/{slavename}.slc_dec',
        slave_par = 'slc_desq/{slavename}.slc_dec.par',
    output:
        off_par = 'diff/{mastername}_{slavename}.off_par',
        ifgram = 'diff/{mastername}_{slavename}.int',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    params:
        mli1 = 'mli/{mastername}.mli',
        mli2 = 'mli/{slavename}.mli'
    run:
        shell("create_offset {input.master_par} {input.slave_par} {output.off_par} - 2 1 0")
        shell("interf_SLC {input.master} {input.slave} {input.master_par} {input.slave_par} {output.off_par} {params.mli1} {params.mli2} {output.ifgram} 2 2 - - 0 0")


rule create_baselines:


rule stacking:
    output:
        avg_ifgram  = 'outputs/{start_dt}_{stop_dt}_{chan}.int'
    input:
        ifgrams = all_ifgrams,
        baselines =  'outputs/{start_dt}_{stop_dt}_{chan}.csv'
        slc_list = 'list_of_slcs.csv'
    wildcard_constraints:
        start_dt = dt_regex,
        stop_dt = dt_regex,
    params:
        ridx = 2502,
        azidx = 60
    run:
        import pyrat.fileutils.gpri_files as gpf
        import pyrat.visualization.visfun as vf
        import numpy as np
        import matplotlib.pyplot as plt
        ifgram = input.ifgrams[0]
        par = ifgram.replace('.int', '.off_par')
        current_if = gpf.gammaDataset(par, ifgram, dtype=gpf.type_mapping['FCOMPLEX'])
        avg_if = np.zeros(current_if.shape,dtype=np.float)
        avg_coh = np.zeros(current_if.shape,dtype=np.float)
        for ifgram in input.ifgrams:
            par = ifgram.replace('.int', '.off_par')
            current_if = gpf.gammaDataset(par, ifgram, dtype=gpf.type_mapping['FCOMPLEX'])
            avg_if += np.angle(current_if * current_if[params.ridx, params.azidx].conj())
            avg_coh += np.abs(current_if)
        avg_coh /= len(input.ifgrams)
        avg_if /= len(input.ifgrams)
        avg_ifgram = avg_coh * np.exp(1j * avg_if)
        rgb, norm, rest = vf.dismph(avg_ifgram, coherence=False)
        plt.imshow(rgb, aspect=1/10)
        plt.show()
