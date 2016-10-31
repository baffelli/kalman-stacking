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

def select_n_dates(string_dates, date_start, n_dates):
    dates = [ dt.datetime.strptime(x, dtfmt) for x in string_dates]
    dt_start = dt.datetime.strptime(date_start, dtfmt)
#    dt_end = dt.datetime.strptime(date_end, dtfmt)
    valid_dates = [date.strftime(dtfmt) for date in dates if dt_start < date][1:int(n_dates)]
    return valid_dates


class StackHelper:
    def __init__(self, list_of_slcs, step=1, stride=1, window=1):
#    """
#        step : int
#            step between master and slave slc
#        stride : int
#            increment of master index
#        window: int
#            the slave interferogram counter is incremented between n_master + 1 and n_master + 1 + window with step size stride
#    """
        with open('list_of_slcs.csv', 'r') as f:#the object contains a list of valid slcs
            reader = csv.reader(f)
            self.all_dates = list(reader)[0]
        self.step = step
        self.stride = stride
        self.window = window


    def itab(self, n_slc):
        tab = []
        for image_counter, idx_master in enumerate(range(1, n_slc, self.step)):
            for idx_slave in range(idx_master + 1, idx_master+1+self.window, self.stride):
                if idx_slave < n_slc:
                    tab.append([idx_master, idx_slave, image_counter, 1 ])
        return tab

#    def slc_tab(self, )


    def valid_dates(self, wildcards):
        return select_date_range(self.all_dates, wildcards.start_dt, wildcards.stop_dt)

    def valid_dates_n(self, wildcards):
        return select_n_dates(self.all_dates, wildcards.start_dt, wildcards.nifgrams)

    def all_diff(self, wildcards):
        valid_dates = self.valid_dates_n(wildcards)
        itab = self.itab(len(valid_dates))
        diffs = []
        for idx_master, idx_slave, *rest in itab:
            current = 'diff/{master}_{chan}_{slave}_{chan}.diff'.format(master=valid_dates[idx_master], slave=valid_dates[idx_slave], chan=wildcards.chan, type=type)
            diffs.append(current)
        return diffs


    def all_mli(self, wildcards):
        return self.all_type(wildcards, 'mli')


    def create_itab(self, output_name):
        with open(output_name, 'w+') as of:
            for line in self.itab:
                of.write(map(str, line))


#    def diff(wildcards)





#    def all_combinations(self, wildcards, ext):
#        #First, find all valid dates
#        valid_dates = self.valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
#        ifgrams = []
#        for vd_master in valid_dates:
#
#        ifgrams = expand(expand('diff/{{vd1}}_{chan}_{{vd2}}_{chan}.{{ext}}', chan=wildcards.chan),zip, vd1=valid_dates, vd2=valid_dates[1::])
#        return ifgrams

#def all_aps(wildcards):
#    valid_dates = valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
#    ifgrams = expand(expand('diff/{{vd1}}_{chan}_{{vd2}}_{chan}.aps', chan=wildcards.chan),zip, vd1=valid_dates, vd2=valid_dates[1::])
#    return ifgrams
#
#
#def plist(wildcards):
#    return expand('ipt/{start_dt}_{stop_dt}_{chan}.plist',start_dt=wildcards.start_dt, stop_dt=wildcards.stop_dt, chan=wildcards.chan)
#
#def msr(wildcards):
#    return expand('ipt/{start_dt}_{stop_dt}_{chan}.msr',start_dt=wildcards.start_dt, stop_dt=wildcards.stop_dt, chan=wildcards.chan)
#
#
#def all_slcs(wildcards):
#    valid_dates = valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
#    slcs = expand('slc_corr/{vd1}_{chan}.slc_dec', chan=wildcards.chan, vd1=valid_dates)
#    return slcs

#def all_mli(wildcards):
#    valid_dates = valid_dates_on_server(wildcards.start_dt, wildcards.stop_dt)
#    mlis = expand('mli/{vd1}_{chan}.mli', chan=wildcards.chan, vd1=valid_dates)
#    return mlis

#
#def wildcards_formatter(input_string, wildcards):
#    wildcard_names = re.match(wildcard_re, str)
#    input_string_no_dot = re.sub(wildcard_re, "\<g>wildcard_name")
#    print(input_string_no_dot)





include: pyrat.rules['raw_to_slc']
include: pyrat.rules['geocoding']


#Initialize stack
stack = StackHelper('list_of_slcs.csv')

rule all:
    input:
        expand('stack/20150803_060519_n_50_AAAl.diff',ext=['diff','int']),
        'list_of_slcs.csv'




##############################################################################
## Correct squint
ruleorder: correct_squint_in_slc > range_compression
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



#Convert shape file to tif containing the area to mask
rule shapefile_to_raster:
    input:
        mask = 'geo/swissTLM3D-2016-tlm_gelaendename_Clip',
        dem = 'geo/swissALTI3D_2016_Clip.tif'
    output:
        mask = 'geo/bisgletscher_mask.tif'
    script:
        'scripts/shapefile_to_raster.py'


#Segment the mask so that it only covers
#the extent covered by DEM
rule segment_mask:
    input:
        mask = 'geo/bisgletscher_mask.tif',
        dem_seg_par = 'geo/' + config['geocoding']['table_name'] + '.dem_seg.par'
    output:
        mask = 'geo/bisgletscher_mask_seg.tif',
        mask_bmp = 'geo/bisgletscher_mask_seg.bmp'
    script:
        'scripts/segment_mask.py'

#Produce a mask of the glacier
rule glacier_mask:
    input:
        mli_ref_par = 'geo/' + config['geocoding']['table_name'] + '.dem_seg.par',
        mask = 'geo/bisgletscher_mask_seg.bmp_fgc'
    output:
        mask = 'geo/bisgletscher_mask.bmp'
    run:
        shell('cp {input.mask} {output.mask}')


#Compute position of reference point
rule reference_coordinate:
    input:
        conf = 'bisgletscher.json',
        dem_par = 'geo/' + config['geocoding']['table_name'] + '.dem_seg.par',
        lut = 'geo/' + config['geocoding']['table_name'] + '.gpri_to_dem',
    output:
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '_reference_pos.csv',
    run:
        import pyrat.geo.geofun as geo
        import csv
        table = geo.GeocodingTable(input.dem_par, input.lut)
        radar_coord = table.geo_coord_to_radar_coord(config['interferogram']['reference_coordinate'])
        with open(output.reference_coord, 'w+') as ouf:
            writer = csv.writer(ouf)
            writer.writerow(radar_coord)



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

##Write slc_tab file
#rule tab:
#    output:
#        'tab/{start_dt}_{stop_dt}_{chan}.slc_tab',
#    input:
#        all_mli
#    wildcard_constraints:
#        start_dt = dt_regex,
#        stop_dt = dt_regex,
#    run:
#        with open(output[0],'w') as ptab:
#            for slc in input:
#                ptab.write(slc + '\t' + slc + '.par\n')

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




ruleorder: mli > multi_look
rule mli:
    input:
        slc = 'slc_desq/{slcname}.slc_dec',
        slc_par = 'slc_desq/{slcname}.slc_dec.par'
    output:
        mli = 'mli/{slcname}.mli',
        mli_par = 'mli/{slcname}.mli.par'
    params:
        rlks = config['interferogram']['rlks'],
        azlks = config['interferogram']['azlks'],
    run:
        shell('multi_look {input.slc} {input.slc_par} {output.mli} {output.mli_par} {params.rlks} {params.azlks} - -')

rule cc:
    input:
        ifgram = 'diff/{mastername}_{slavename}.int',
        ifgram_par = 'diff/{mastername}_{slavename}.int_par',
        mli1  = 'mli/{mastername}.mli',
        mli2  = 'mli/{slavename}.mli',
        mli1_par  = 'mli/{mastername}.mli.par',
    output:
        cc = 'diff/{mastername}_{slavename}.cc',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
       import pyrat.fileutils.gpri_files as gpf
       wd = gpf.par_to_dict(input.mli1_par)['range_samples']
       cc_cmd = "cc_wave {{input.ifgram}} {{input.mli1}} {{input.mli2}} {{output.cc}} {wd} - - 0".format(wd=wd)
       shell(cc_cmd)

rule cc_mask:
    output:
        cc_mask = 'diff/{mastername}_{slavename}.cc_mask.bmp',
    input:
        cc =  'diff/{mastername}_{slavename}.cc',
        pwr = 'mli/{mastername}.mli',
        mli_par = 'mli/{mastername}.mli.par',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
       import pyrat.fileutils.gpri_files as gpf
       wd = gpf.par_to_dict(input.mli_par)['range_samples']
       mask_cmd = "rascc_mask {{input.cc}} {{input.pwr}} {wd} - - - - - 0.5 0.01 - - - - - {{output.cc_mask}}".format(wd=wd)
       shell(mask_cmd)

#Invert a mask
rule invert_mask:
    input:
        mask = "{name}_mask.bmp"
    output:
        mask_inv = "{name}_mask_inv.bmp"
    run:
        print(wildcards)
        shell('mask_op {input.mask} {input.mask} {output.mask_inv} 2')

#Convert a layover/shadow map into a bitmap
rule ls_map_bmp:
    input:
        mask = 'geo/' + config['geocoding']['table_name'] + '.sh_map_fgc',
        ref_mli_par = config['geocoding']['ref_mli_par']
    output:
        bmp_mask = "geo/bisgletscher_shadow_layover.bmp"
    script:
        'scripts/layover_shadow_as_bmp.py'

rule glacier_validity_mask:
    input:
        cc_mask = 'diff/{mastername}_{slavename}.cc_mask.bmp',
        glacier_mask = 'geo/bisgletscher_mask.bmp',
    output:
        validity_mask = 'diff/{mastername}_{slavename}.unw_mask.bmp'
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
        shell("mask_op {input.cc_mask} {input.glacier_mask} {output.validity_mask} 0")

rule ifgram:
    input:
        master = 'slc_desq/{mastername}.slc_dec',
        master_par = 'slc_desq/{mastername}.slc_dec.par',
        slave = 'slc_desq/{slavename}.slc_dec',
        slave_par = 'slc_desq/{slavename}.slc_dec.par',
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '_reference_pos.csv',
    output:
        int_par = 'diff/{mastername}_{slavename}.int_par',
        ifgram = temp('diff/{mastername}_{slavename}.int_unref'),
        ifgram_ref = 'diff/{mastername}_{slavename}.int',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    params:
        mli1 = 'mli/{mastername}.mli',
        mli2 = 'mli/{slavename}.mli',
        rlks = config['interferogram']['rlks'],
        azlks = config['interferogram']['azlks'],
    run:
        import numpy as np
        import pyrat.fileutils.gpri_files as gpf
        import snakemake
        shell("create_offset {input.master_par} {input.slave_par} {output.int_par} - 2 1 0")
        shell(
        "SLC_intf {input.master} {input.slave} {input.master_par} {input.slave_par} {output.int_par} {output.ifgram} {params.rlks} {params.azlks} - - 0 0 1 - - - - -")
        ifgram_par = gpf.par_to_dict(output.int_par)
        # Load reference coord
        ref_coord = np.genfromtxt(input.reference_coord, delimiter=',')
        print(ref_coord)
        # referencing
        ref_cmd = "cpx_math {{output.ifgram}} - {{output.ifgram_ref}} {wd} 0 {ridx} {azidx} {nr} {naz} - - - 1".format(
        ridx=ref_coord[0], azidx=ref_coord[1], nr=config['interferogram']['reference_region_size'][0],
        naz=config['interferogram']['reference_region_size'][1], wd=ifgram_par.interferogram_width)
        shell(ref_cmd)

#Compute the atmospheric phase screen for an image
rule aps:
    output:
        aps_unw= 'diff/{mastername}_{slavename}.aps_unw',
        aps = 'diff/{mastername}_{slavename}.aps',
        int_filt = temp('diff/{mastername}_{slavename}.int_filt'),
        int_mask = temp('diff/{mastername}_{slavename}.int_masked'),
    input:
        ifgram = 'diff/{mastername}_{slavename}.int',
        int_par = 'diff/{mastername}_{slavename}.int_par',
        mask = "diff/{mastername}_{slavename}.unw_mask.bmp",
        mli_par = 'mli/{mastername}.mli.par',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    params:
        aps_window = 100,
        filter_type = 1
    run:
        import pyrat.fileutils.gpri_files as gpf
        wd = gpf.par_to_dict(input.mli_par)['range_samples']
        mask_cmd = "mask_data {{input.ifgram}} {wd} {{output.int_mask}} {{input.mask}} 1".format(wd=wd)
        shell(mask_cmd)
        filt_cmd = "fspf {{output.int_mask}} {{output.int_filt}} {wd} 0 {{params.aps_window}} 3 {{input.mli_par}}".format(wd=wd)
        shell(filt_cmd)
        mcf_cmd = "mcf {{output.int_filt}} - - {{output.aps_unw}} {wd} - - - - - 1 1 - - - 0 ".format(wd=wd)
        shell(mcf_cmd)
        #interpolate
        interp_cmd = "interp_ad {{output.aps_unw}} {{output.aps}} {wd} 250 10 250 3 2".format(wd=wd)
        shell(interp_cmd)


rule cleanup_diff:
    run:
        shell('rm diff/*.aps')
        shell('rm diff/*.diff*')
        shell('rm diff/*cc_mask*')
        shell('rm diff/*unw_mask*')

rule cleanup_geo:
    run:
        shell('rm geo/Dom*')
        shell('rm geo/bisgletscher*')


#Remove the aps from the interferogram
rule diff_ifgram:
    output:
        diff_int = 'diff/{mastername}_{slavename}.diff',
        diff_par = 'diff/{mastername}_{slavename}.diff_par',
    input:
        aps = 'diff/{mastername}_{slavename}.aps',
        int = 'diff/{mastername}_{slavename}.int',
        int_par = 'diff/{mastername}_{slavename}.int_par',
        mli1_par = "mli/{mastername}.mli.par",
        mli2_par=  "mli/{slavename}.mli.par",
        mask = "diff/{mastername}_{slavename}.unw_mask.bmp",
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
#        #Create diff_par:
        par_cmd = "create_diff_par {input.int_par} - {output.diff_par} 0 0"
        shell(par_cmd)
        from pyrat.diff.core import Interferogram as intgram
        import pyrat.fileutils.gpri_files as gpf
        import numpy as np
        ifgram = intgram(output.diff_par, input.int, master_par = input.mli1_par, slave_par = input.mli2_par, dtype=gpf.type_mapping['FCOMPLEX'])
        aps = gpf.gammaDataset(output.diff_par, input.aps, dtype=gpf.type_mapping['FLOAT'])
        ifgram = np.exp(-1j * np.array(aps)) * ifgram
        ifgram.tofile(output.diff_par, output.diff_int)


#Projects a single corrected interferogram to a map
rule diff_to_map:
    output:
        map = 'outputs/{filename}.diff.pdf'
    input:
        diff = 'diff/{filename}.diff_gc.tif',
        basemap = 'geo/pk25krel_latest_Clip.tif',
        mask = 'geo/swissTLM3D-2016-tlm_gelaendename_Clip',
    script:
        'scripts/overlay_displacement_on_map.py'





rule stacking:
    output:
        avg_ifgram  = 'stack/{start_dt}_n_{nifgrams}_{chan}.diff',
        avg_ifgram_par  = 'stack/{start_dt}_n_{nifgrams}_{chan}.diff_par',
    input:
        ifgrams = stack.all_diff,
#        mli = stack.all_mli,
        slc_list = 'list_of_slcs.csv',
        lut = 'geo/' + config['geocoding']['table_name']  + '.gpri_to_dem',
        dem_par = 'geo/' + config['geocoding']['table_name']  + '.dem_seg.par',
        map = 'geo/pk25krel_latest_Clip_24bt.tif_fgc'
    wildcard_constraints:
        start_dt = dt_regex,
#        stop_dt = dt_regex,
    params:
        ridx = 1173,
        azidx = 112,
        ws = 8
    script:
        'scripts/stacking.py'


