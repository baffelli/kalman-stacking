import pyrat
from snakemake.remote.SFTP import RemoteProvider, RemoteObject
import datetime as dt
import csv
import re
provider = RemoteProvider(port=22,username="baffelli", private_key="/home/baffelli/.ssh/id_rsa")

#tar = provider.glob_wildcards("ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/raw/{date,20150803}_{time,06\d+}.tar")

configfile: './bisgletscher.json'
slc_regex = "(\d{8})_(\d{6})((_[A-B]{3}[l-u]))?"
dt_regex = "(\d{8})_(\d{6})"
dtfmt = "%Y%m%d_%H%M%S"
wildcard_re = "\{\S+\.(?<wildcard_name>\S+)\}"

def select_date_range(string_dates, date_start, date_end):
    dates = [ dt.datetime.strptime(x, dtfmt) for x in string_dates]
    dt_start = dt.datetime.strptime(date_start, dtfmt)
    dt_end = dt.datetime.strptime(date_end, dtfmt)
    valid_dates = [date.strftime(dtfmt) for date in dates if dt_start < date < dt_end]
    return valid_dates

def valid_dates_on_server(start_dt, stop_dt):
    with open('list_of_slcs.csv', 'r') as f:
        reader = csv.reader(f)
        all_dates = list(reader)
    valid_dates = select_date_range(all_dates[0], start_dt, stop_dt)
    return valid_dates

def all_ifgrams(wildcards):
    valid_dates = valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
    ifgrams = expand(expand('diff/{{vd1}}_{chan}_{{vd2}}_{chan}.int', chan=wildcards.chan),zip, vd1=valid_dates, vd2=valid_dates[1::])
    return ifgrams

def plist(wildcards):
    return expand('ipt/{start_dt}_{stop_dt}_{chan}.plist',start_dt=wildcards.start_dt, stop_dt=wildcards.stop_dt, chan=wildcards.chan)


def all_slcs(wildcards):
    valid_dates = valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
    slcs = expand('slc_corr/{vd1}_{chan}.slc_dec', chan=wildcards.chan, vd1=valid_dates)
    return slcs

def all_mli(wildcards):
    valid_dates = valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
    mlis = expand('mli/{vd1}_{chan}.mli', chan=wildcards.chan, vd1=valid_dates)
    return mlis


def wildcards_formatter(input_string, wildcards):
    wildcard_names = re.match(wildcard_re, str)
    input_string_no_dot = re.sub(wildcard_re, "\<g>wildcard_name")
    print(input_string_no_dot)





include: pyrat.rules['raw_to_slc']


rule all:
    input:
#       'outputs/20150803_060019_20150803_063519_AAAl_BBBl.int.pdf',
       'outputs/20150803_060019_20150803_063519_AAAl_BBBl.int',
       'list_of_slcs.csv'




##############################################################################
## Correct squint
rule correct_squint_in_slc:
		input:
			slc = "slc_chan/{dataname}_{chan,[A-B]{3}}{rx}.slc",
			slc_par = "slc_chan/{dataname}_{chan,[A-B]{3}}{rx}.slc.par",
		output:
			corr =  "slc_desq/{dataname}_{chan,[A-B]{3}}{rx}.slc",
			corr_par = "slc_desq/{dataname}_{chan,[A-B]{3}}{rx}.slc.par"
		params:
			squint_rate = lambda wildcards: config["desquint"]["{chan}_squint_rate".format(chan=wildcards.chan)],
		run:
			slc = gpf.gammaDataset(input.slc_par, input.slc)
			raw_desq = gpf.correct_squint_in_SLC(slc, squint_rate=params.squint_rate)
			slc.tofile(output.corr_par, output.corr)


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

#Write slc_tab file
rule tab:
    output:
        'tab/{start_dt}_{stop_dt}_{chan}.slc_tab',
    input:
        all_mli
    wildcard_constraints:
        start_dt = dt_regex,
        stop_dt = dt_regex,
    run:
        with open(output[0],'w') as ptab:
            for slc in input:
                ptab.write(slc + '\t' + slc + '.par\n')

#Average power
rule ave_pwr:
    output:
        ave ='mli/{start_dt}_{stop_dt}_{chan}.ave_pwr'
    input:
        tab = 'tab/{start_dt}_{stop_dt}_{chan}.slc_tab',
    wildcard_constraints:
        start_dt = dt_regex,
        stop_dt = dt_regex,
    run:
        import pyrat.fileutils.gpri_files as gpf
        with open(input.tab) as tab:
            mli, mli_par = tab.readline().replace('\n','').split('\t')
            print(mli_par)
        width = gpf.get_width(mli_par)
        ave_cmd = "ave_image {{input.tab}} {width} {{output.ave}} - - - - 0".format(width=width)
        shell(ave_cmd)


#create candidate point target list
rule ptarg_screen:
    output:
         msr = 'ipt/{start_dt}_{stop_dt}_{chan}.msr',
         plist = 'ipt/{start_dt}_{stop_dt}_{chan}.plist',
    wildcard_constraints:
        start_dt = dt_regex,
        stop_dt = dt_regex,
    input:
         tab = 'tab/{start_dt}_{stop_dt}_{chan}.slc_tab',
         ave_pwr = 'mli/{start_dt}_{stop_dt}_{chan}.ave_pwr'
    run:
        with open(input.tab) as tab:
            slc, slc_par = tab.readline().replace('\n','').split('\t')
        stat_cmd = "pwr_stat {{input.tab}} {slc_par} {{output.msr}} {{output.plist}} - - - - - -".format(slc_par = slc_par)
        shell(stat_cmd)

ruleorder: mli > multi_look
rule mli:
    input:
        slc = 'slc_corr/{slcname}.slc_dec',
        slc_par = 'slc_corr/{slcname}.slc_dec.par'
    output:
        mli = 'mli/{slcname}.mli',
        mli_par = 'mli/{slcname}.mli.par'
    params:
        rlks = config['interferogram']['rlks'],
        azlks = config['interferogram']['azlks'],
    run:
        shell('multi_look {input.slc} {input.slc_par} {output.mli} {output.mli_par} {params.rlks} {params.azlks}')

rule cc:
    input:
        ifgram = 'diff/{mastername}_{slavename}.int',
        ifgram_par = 'diff/{mastername}_{slavename}.int_par',
        mli1  = 'mli/{mastername}.mli',
        mli2  = 'mli/{slavename}.mli',
    output:
        cc = 'diff/{mastername}_{slavename}.int',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
       import pyrat.fileuitls.gpri_files as gpf




rule ifgram:
    input:
        master = 'slc_corr/{mastername}.slc_dec',
        master_par = 'slc_corr/{mastername}.slc_dec.par',
        slave = 'slc_corr/{slavename}.slc_dec',
        slave_par = 'slc_corr/{slavename}.slc_dec.par',
    output:
        int_par = 'diff/{mastername}_{slavename}.int_par',
        ifgram = 'diff/{mastername}_{slavename}.int',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    params:
        mli1 = 'mli/{mastername}.mli',
        mli2 = 'mli/{slavename}.mli',
        rlks = config['interferogram']['rlks'],
        azlks = config['interferogram']['azlks'],
    run:
        shell("create_offset {input.master_par} {input.slave_par} {output.int_par} - 2 1 0")
        shell("SLC_intf {input.master} {input.slave} {input.master_par} {input.slave_par} {output.int_par} {output.ifgram} {params.rlks} {params.azlks} - - 0 0 1 - - - - -")
        #Compute temporal baseline
        par1 = gpf.par_to_dict(input.master_par)
        par2 = gpf.par_to_dict(input.slave_par)
        int_par = gpf.par_to_dict(output.int_par)
        bl = par2['center_time'][0] - par1['center_time'][0]
        int_par['temporal_baseline'] = [bl, 's']
        gpf.dict_to_par(int_par, output.int_par)





rule create_baselines:
    output:
        baselines =  'outputs/{start_dt}_{stop_dt}_{chan}.csv'
    input:
        ifgrams = all_ifgrams,
    run:
        import re
        import csv
        par_str = "slc_corr/{slcname}.slc.par"
        with open(output.baselines) as of:
            writer= csv.writer(of)
            for ifgram in ifgrams:
                slc_1, slc_2 = re.split(slc_regex, ifgram)
                par1 = gpf.par_to_dict(str.format(slc_1))
                par2 = gpf.par_to_dict(par_str.format(slc_2))
                bl = par1['center_time'][1] - par1['center_time'][0]
                out_str = [slc_1, slc_2, bl]
                write.writerow(out_str)


rule stacking:
    output:
        avg_ifgram  = 'stack/{start_dt}_{stop_dt}_{chan}.int',
        avg_ifgram_par  = 'stack/{start_dt}_{stop_dt}_{chan}.int_par',
    input:
        ptarg_screen = plist,
        ifgrams = all_ifgrams,
        slc_list = 'list_of_slcs.csv'
    wildcard_constraints:
        start_dt = dt_regex,
        stop_dt = dt_regex,
    params:
        ridx = 1675,
        azidx = 123,
        ws = 8
    run:
        import pyrat.fileutils.gpri_files as gpf
        import pyrat.visualization.visfun as vf
        import numpy as np
        import matplotlib.pyplot as plt
        import pyrat.core.corefun
        ifgram = input.ifgrams[0]
        par = ifgram + "_par"
        current_if = gpf.gammaDataset(par, ifgram, dtype=gpf.type_mapping['FCOMPLEX'])
        avg_if = np.zeros(current_if.shape,dtype=np.float)
        avg_coh = np.zeros(current_if.shape,dtype=np.float)
        dt_avg = 0
        for ifgram in input.ifgrams:
            par = ifgram + '_par'
            current_if = gpf.gammaDataset(par, ifgram, dtype=gpf.type_mapping['FCOMPLEX'])
#            current_if = cf.smooth(current_if, [3,1])
            dt = current_if.temporal_baseline[0]
            #extract phase of refernce region
            ref_region = cf.window_idx(current_if, [params.ridx, params.azidx], [params.ws, params.ws])
            ref_phase = np.mean(np.angle(current_if[ref_region]))
            avg_if += np.angle(current_if * ref_phase.conj()) * dt
            avg_coh += np.abs(current_if) * dt
            dt_avg += dt
        avg_coh /= dt_avg
        avg_if /= dt_avg
        avg_ifgram = avg_coh * np.exp(1j * avg_if)
        #TODO fix write_dataset so that it ssaves the correct type
        gpf.write_dataset(avg_ifgram.astype(gpf.type_mapping['FCOMPLEX']), current_if.__dict__, output.avg_ifgram_par, output.avg_ifgram)

rule compare_pol:
    input:
        avg_ifgram_1  = 'stack/{start_dt}_{stop_dt}_{chan1}.int',
        avg_ifgram_2  = 'stack/{start_dt}_{stop_dt}_{chan2}.int',
        avg_ifgram_1_par  = 'stack/{start_dt}_{stop_dt}_{chan1}.int_par',
        avg_ifgram_2_par  = 'stack/{start_dt}_{stop_dt}_{chan2}.int_par',
    output:
        diff = 'outputs/{start_dt}_{stop_dt}_{chan1}_{chan2}.int',
    wildcard_constraints:
        start_dt = dt_regex,
        stop_dt = dt_regex
    run:
        import pyrat.fileutils.gpri_files as gpf
        import pyrat.visualization.visfun as vf
        import numpy as np
        import matplotlib.pyplot as plt
        HH = gpf.gammaDataset(input.avg_ifgram_1_par, input.avg_ifgram_1, dtype=gpf.type_mapping['FCOMPLEX'])
        VV = gpf.gammaDataset(input.avg_ifgram_2_par, input.avg_ifgram_2, dtype=gpf.type_mapping['FCOMPLEX'])
        rgb_HH, norm, rest = vf.dismph(HH, sf = 1, k=0.7)
        rgb_VV, norm, rest = vf.dismph(VV, sf = 1, k=0.7)
        ax = plt.subplot(1,3,1)
        plt.imshow(rgb_HH, aspect=1/10)
        plt.subplot(1,3,2, sharex=ax, sharey=ax)
        plt.imshow(rgb_VV,  aspect=1/10)
        plt.subplot(1,3,3, sharex=ax, sharey=ax)
        plt.imshow(np.angle(VV * HH.conj()), aspect=1/10, cmap='jet')
        plt.show()
