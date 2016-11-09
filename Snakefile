import pyrat
from snakemake.remote.SFTP import RemoteProvider, RemoteObject
import datetime as dt
import csv
import re
provider = RemoteProvider(port=22,username="baffelli", private_key="/home/baffelli/.ssh/id_rsa")



configfile: './bisgletscher.json'
slc_regex = r"(?P<dt>\d{8}_\d{6})(_[A-B]{3}[l-u])?"
dt_regex = "(\d{8})_(\d{6})"
dtfmt = "%Y%m%d_%H%M%S"
wildcard_re = "\{\S+\.(?<wildcard_name>\S+)\}"


def select_date_range(string_dates, date_start, date_end):
    dates = [ dt.datetime.strptime(x, dtfmt) for x in string_dates]
    dt_start = dt.datetime.strptime(date_start, dtfmt)
    dt_end = dt.datetime.strptime(date_end, dtfmt)
    valid_dates = [date.strftime(dtfmt) for date in dates if dt_start < date < dt_end]
    return valid_dates

def select_nth_previous(string_dates, date_end, n):
    dates = [ dt.datetime.strptime(x, dtfmt) for x in string_dates]
    dt_end = dt.datetime.strptime(date_end, dtfmt)
    #find the closest element to the end date
    closest_index = min(dates, key=lambda x: abs(x- dt_end)) - 1
    start = closest_index - n
    valid_dates = dates[start]
    return valid_dates

def select_n_dates(string_dates, date_start, n_dates):
    dates = [ dt.datetime.strptime(x, dtfmt) for x in string_dates]
    dt_start = dt.datetime.strptime(date_start, dtfmt)
#    dt_end = dt.datetime.strptime(date_end, dtfmt)
    valid_dates = [date.strftime(dtfmt) for date in dates if dt_start < date][1:int(n_dates)]
    return valid_dates

def get_reference_coord(dict):
    reference_feature = [f for f in dict['features'] if f['id']=='reference'][0]
    radar_coord = reference_feature['properties']['radar_coordinates']
    return radar_coord

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
        with open(list_of_slcs, 'r') as f:#the object contains a list of valid slcs
            reader = csv.reader(f)
            self.all_dates = list(reader)[0]
        self.step = step
        self.stride = stride
        self.window = window


    def itab(self, n_slc, window, stride, step):
        tab = []
        for image_counter, idx_master in enumerate(range(1, n_slc, step)):
            for idx_slave in range(idx_master + 1, idx_master+1+window, stride):
                if idx_slave < n_slc:
                    tab.append([idx_master, idx_slave, image_counter, 1 ])
        return tab


    def valid_dates(self, wildcards):
        return select_date_range(self.all_dates, wildcards.start_dt, wildcards.stop_dt)

    def valid_dates_n(self, wildcards):
        return select_n_dates(self.all_dates, wildcards.start_dt, wildcards.nifgrams)

    def previous_stack(self, wildcards):
        """
            Return the dates of the n slcs preceding the current one
        """
        matches = re.search(slc_regex, wildcards.slavename)
        stop_dt = matches.group(0)
        start =  select_nth_previous(self.all_dates, stop_dt, 20)
        return 'stack/{start}_stacking_20.diff'.format(start=start)

    def all_pairs(self, pattern, wildcards):
        valid_dates = self.valid_dates_n(wildcards)
        itab = self.itab(len(valid_dates))
        diffs = []
        for idx_master, idx_slave, *rest in itab:
            current = pattern.format(master=valid_dates[idx_master], slave=valid_dates[idx_slave], chan=wildcards.chan, type=type)
            diffs.append(current)
        return diffs

    def all_diff(self, wildcards):
        return self.all_pattern('diff/{master}_{chan}_{slave}_{chan}.diff', wildcards)

    def first_valid_mli_par(self, wildcards):
        #"Return the first valid mli par from a list given the wildcards"
        valid_dates = self.valid_dates_n(wildcards)
        return 'mli/{date}_{chan}.mli.par'.format(date=valid_dates[0], chan=wildcards.chan)



    def all_single(self, pattern, wildcards):
        valid_dates = self.valid_dates_n(wildcards)
        files = expand(pattern, date=valid_dates, **wildcards)
        return files

    def all_slc(self, wildcards):
        valid_dates = self.valid_dates_n(wildcards)
        files = expand('slc_desq/{date}_{chan}.slc_dec', date=valid_dates, chan=wildcards.chan)

    def stacking_inputs(self, wildcards):
        pref_mapping = {'mli':'mli', 'cc':'int', 'diff':'diff', 'unw':'diff', 'aps':'diff'}
        function_map = {'mli':self.all_single, 'int':self.all_pairs, 'diff':self.all_pairs}
        pattern_map = {'mli':'mli/{date}_{chan}.{ext}', 'int':"int/{{master}}_{{chan}}_{{slave}}_{{chan}}.{ext}", 'diff':"diff/{{master}}_{{chan}}_{{slave}}_{{chan}}.{ext}"}
        pref = pref_mapping[wildcards.ext]
        function = function_map[pref]
        pattern = pattern_map[pref].format(ext=wildcards.ext)
        inputs = function(pattern, wildcards)
        return inputs


    def create_itab(self, n_slc,output_name, ):
        with open(output_name, 'w+') as of:
            itab=self.itab(n_slc)
            for line in itab:
                of.write(map(str, line))








include: pyrat.rules['raw_to_slc']
include: pyrat.rules['geocoding']


#Initialize stack
stack = StackHelper('list_of_slcs.csv')

#Find all ifgrams
ints, = glob_wildcards('int/{name}.int')
ras = expand('int/{ints}.int.bmp', ints=ints)

rule all:
    input:
#        expand('stack/20150803_060519_stack_{n}_AAAl.{ext}',ext=['cc.ave_gc.tif', 'unw.ave_gc.tif', 'diff.ave_gc.tif'], n=[10,20]),
#        expand("diff/20150803_120249_AAAl_20150803_120519_AAAl.{ft}",ft=['aps_ref']),
#        expand('stack/20150803_060519_stacl_{n}_AAAl.variogram', n=[10,20]),
        expand('ipt/20150803_060519_stack_20_AAAl.{ext}', ext=['paps', 'plist_glacier', 'phgt_masked', 'phgt_glacier']),
        "mli/20150803_060749_AAAl.mli",
        'geo/Dom.ls_map.tif'

#        'list_of_slcs.csv'



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
			print(type(raw_desq))
			raw_desq.tofile(output.corr_par, output.corr)



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
        3
        shell('cp {input.mask} {output.mask}')


#Compute position of reference point and store it in json file
rule reference_coordinate:
    input:
        dem_par = 'geo/' + config['geocoding']['table_name'] + '.dem_seg.par',
        lut = 'geo/' + config['geocoding']['table_name'] + '.gpri_to_dem',
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '.shp',
    output:
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '.json',
    script:
        'scripts/shapefile_to_json.py'


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
        print('doing')
        slcs = provider.glob_wildcards("ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/raw/{slcname}.tar")
        with open(output[0], 'w+') as of:
            writer = csv.writer(of)
            writer.writerows(slcs)

ruleorder: multi_look > average
#Average unwrapped data/coherence/aps etc
rule average:
    output:
        ave ='stack/{start_dt}_stack_{nifgrams}_{chan}.{ext}.ave',
        tab = 'stack/{start_dt}_stack_{nifgrams}_{chan}.{ext}.tab'
    input:
        tab = stack.stacking_inputs,
        mli_par = stack.first_valid_mli_par,
    wildcard_constraints:
        start_dt = dt_regex,
    run:
        import pyrat.fileutils.gpri_files as gpf
        with open(output.tab, 'w+') as tab:
            for file in input.tab:
                print(file)
                tab.write(str(file) + '\n')
            print(tab.readlines())
        width = gpf.get_width(input.mli_par)
        ave_cmd = "ave_image {{output.tab}} {width} {{output.ave}} - - - - 1".format(width=width)
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
        ifgram = 'int/{mastername}_{slavename}.int',
        ifgram_par = 'int/{mastername}_{slavename}.int_par',
        mli1  = 'mli/{mastername}.mli',
        mli2  = 'mli/{slavename}.mli',
        mli1_par  = 'mli/{mastername}.mli.par',
    output:
        cc = 'int/{mastername}_{slavename}.cc',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    params:
        rlks = config['interferogram']['rlks'],
        azlks = config['interferogram']['azlks']
    run:
       import pyrat.fileutils.gpri_files as gpf
       wd = gpf.par_to_dict(input.mli1_par)['range_samples']
       cc_cmd = "cc_wave {{input.ifgram}} {{input.mli1}} {{input.mli2}} {{output.cc}} {wd} {{params.rlks}} {{params.azlks}} 0".format(wd=wd)
       shell(cc_cmd)

rule cc_mask:
    output:
        cc_mask = 'diff/{mastername}_{slavename}.cc_mask.bmp',
    input:
        cc =  'int/{mastername}_{slavename}.cc.sm',
        pwr = 'mli/{mastername}.mli',
        mli_par = 'mli/{mastername}.mli.par',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    params:
        coherence_threshold = config['unwrapping']['coherence_threshold']
    run:
       import pyrat.fileutils.gpri_files as gpf
       wd = gpf.par_to_dict(input.mli_par)['range_samples']
       mask_cmd = "rascc_mask {{input.cc}} {{input.pwr}} {wd} - - - - - {{params.coherence_threshold}} - - - - - - {{output.cc_mask}}".format(wd=wd)
       shell(mask_cmd)

#Invert a mask
rule invert_mask:
    input:
        mask = "{name}.bmp"
    output:
        mask_inv = "{name}_inv.bmp"
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
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '.json',
    output:
        int_par = 'int/{mastername}_{slavename}.int_par',
        ifgram = temp('int/{mastername}_{slavename}.int_unref'),
        ifgram_ref = 'int/{mastername}_{slavename}.int',
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
        import json
        shell("create_offset {input.master_par} {input.slave_par} {output.int_par} - {params.rlks} {params.azlks} 0")
        shell(
        "SLC_intf {input.master} {input.slave} {input.master_par} {input.slave_par} {output.int_par} {output.ifgram} {params.rlks} {params.azlks} - - 0 0 1 1 - - - -")
        ifgram_par = gpf.par_to_dict(output.int_par)
        master_par = gpf.par_to_dict(input.master_par)
        slave_par = gpf.par_to_dict(input.slave_par)
        bl = master_par.start_time - slave_par.start_time
        ifgram_par.add_parameter('temporal_baseline', bl, unit='s')
        gpf.dict_to_par(ifgram_par, output.int_par)
        # Load reference coord
        with open(input.reference_coord) as inputfile:
            ref_coord =  get_reference_coord(json.load(inputfile))
        # referencing
        ref_cmd = "cpx_math {{output.ifgram}} - {{output.ifgram_ref}} {wd} 0 {ridx} {azidx} {nr} {naz} - - - 1".format(
        ridx=ref_coord[0], azidx=ref_coord[1], nr=config['interferogram']['reference_region_size'][0],
        naz=config['interferogram']['reference_region_size'][1], wd=ifgram_par.interferogram_width)
        shell(ref_cmd)
        print("Temporal baseline:{bl}".format(bl=bl))


rule rasmph:
    output:
        ras = 'int/{mastername}_{slavename}.int.bmp',
    input:
        master_mli = 'mli/{mastername}.mli',
        int_par = 'int/{mastername}_{slavename}.int_par',
        ifgram = 'int/{mastername}_{slavename}.int.sm',
        cc = 'int/{mastername}_{slavename}.cc.sm',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
        import pyrat.fileutils.gpri_files as gpf
        wd = gpf.par_to_dict(input.int_par)['interferogram_width']
        rascmd = "rasmph_pwr {{input.ifgram}} {{input.master_mli}} {wd} - - - - - - - - {{output.ras}} {{input.cc}}".format(wd=wd)
        shell(rascmd)


#Smooth interferogram
rule adf:
    input:
        ifgram  = 'int/{mastername}_{slavename}.int',
        mli_par = 'mli/{mastername}.mli.par',
    output:
        int_sm = 'int/{mastername}_{slavename}.int.sm',
        cc = 'int/{mastername}_{slavename}.cc.sm',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
        import pyrat.fileutils.gpri_files as gpf
        wd = gpf.par_to_dict(input.mli_par)['range_samples']
        adf_cmd = "adf {{input.ifgram}} {{output.int_sm}} {{output.cc}} {wd} {{config[interferogram][adf_alpha]}} {{config[interferogram][adf_window]}}".format(wd=wd)
        shell(adf_cmd)

##Compute unwrapped interferogram
rule unwrap:
    input:
        ifgram = 'int/{mastername}_{slavename}.int.sm',
        int_par = 'int/{mastername}_{slavename}.int_par',
        mask = "diff/{mastername}_{slavename}.unw_mask.bmp",
        cc_mask = "diff/{mastername}_{slavename}.cc_mask.bmp",
        cc = "int/{mastername}_{slavename}.cc.sm",
        mli_par = 'mli/{mastername}.mli.par',
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '.json',
    output:
        unw = 'diff/{mastername}_{slavename}.unw',
    params:
        mode = config['interferogram']['triangulation_mode']
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
        import pyrat.fileutils.gpri_files as gpf
        import numpy as np
        #read location of reference point
        with open(input.reference_coord) as inputfile:
            ref_coord =  get_reference_coord(json.load(inputfile))
        #get width of data
        wd = gpf.par_to_dict(input.mli_par)['range_samples']
        #unwrap data
        mcf_cmd = "mcf {{input.ifgram}} {{input.cc}} {{input.cc_mask}} {{output.unw}} {wd} {{params.mode}} - - - - 1 1 - {ridx} {azidx} 1 ".format(wd=wd, ridx=ref_coord[0], azidx=ref_coord[1])
        shell(mcf_cmd)


#Create differential interferogram parameters
rule diff_par:
    output:
        diff_par = 'diff/{mastername}_{slavename}.diff_par'
    input:
        int_par = 'int/{mastername}_{slavename}.int_par'
    run:
        par_cmd = "create_diff_par {input.int_par} - {output.diff_par} 0 0"
        shell(par_cmd)

#Compute the atmospheric phase screen for an image
rule aps:
    output:
#        aps_unw= 'diff/{mastername}_{slavename}.aps_unw',
        aps = 'diff/{mastername}_{slavename}.aps',
#        aps_masked = temp('diff/{mastername}_{slavename}.aps.masked'),
##        int_filt = temp('diff/{mastername}_{slavename}.int_filt'
#        diff_par = 'diff/{mastername}_{slavename}.diff_par',
        int_filt = temp('diff/{mastername}_{slavename}.unw.filt'),
        int_masked = temp('diff/{mastername}_{slavename}.unw.mask'),
#        int_unw = temp('diff/{mastername}_{slavename}.int.unw'),
    input:
        ifgram = 'diff/{mastername}_{slavename}.unw',
        int_par = 'int/{mastername}_{slavename}.int_par',
        combined_mask = "diff/{mastername}_{slavename}.unw_mask.bmp",
        cc_mask = "diff/{mastername}_{slavename}.cc_mask.bmp",
        cc = "int/{mastername}_{slavename}.cc",
        mli_par = 'mli/{mastername}.mli.par',
        topo = 'geo/Dom.dem_seg_fgc',
        diff_par = 'diff/{mastername}_{slavename}.diff_par',
        glacier_mask = 'geo/bisgletscher_mask.bmp',
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '.json',
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    params:
        filter_type = 1,
    run:
        import pyrat.fileutils.gpri_files as gpf
        import numpy as np
        wd = gpf.par_to_dict(input.mli_par)['range_samples']
        with open(input.reference_coord) as inputfile:
            ref_coord =  get_reference_coord(json.load(inputfile))
#        #Mask the glacier
        mask_cmd = "mask_data {{input.ifgram}} {wd} {{output.int_masked}} {{input.combined_mask}} 0".format(wd=wd)
        shell(mask_cmd)
        #Filter strongly
        filt_cmd = "fspf {{input.ifgram}} {{output.int_filt}} {wd} 2 {{config[aps][filter_radius]}} 3 {{input.mli_par}}".format(wd=wd)
        shell(filt_cmd)
#        mcf_cmd = "mcf {{output.int_filt}} {{input.cc}} {{input.combined_mask}} {{output.int_unw}} {wd} - - - - - 1 1 - {ridx} {azidx} 1 ".format(wd=wd, ridx=ref_coord[0], azidx=ref_coord[1])
#        shell(mcf_cmd)
        interp_cmd = "interp_ad {{output.int_filt}} {{output.aps}} {wd} 100 10 100 3 2".format(wd=wd)
        shell(interp_cmd)
        #Mask unwrapped
#        mask_cmd = "mask_data {{input.ifgram}} {wd} {{output.int_masked}} {{input.mask}} 1".format(wd=wd)
        #Atmospheric model
#        mod_cmd = "atm_mod {{input.ifgram}} {{input.topo}} {{input.diff_par}} {{output.aps}} - - {{input.glacier_mask}} 1 {ridx} {azidx}".format(ridx=int(ref_coord[0]), azidx=int(ref_coord[1]))
#        shell(mod_cmd)
##        shell(mod_cmd)
#        #low_pass filter
#        filt_cmd = "fspf {{output.int_masked}} {{output.int_interp}} {wd} 2 {{params.aps_window}} 3 {{input.mli_par}}".format(wd=wd)
##        #Close holes
#        interp_cmd = "interp_ad {{output.int_interp}} {{output.aps}} {wd} 100 10 100 3 2".format(wd=wd)
#        shell(mask_cmd)
#        shell(filt_cmd)
#        shell(interp_cmd)
#        #load reference coordinate
#        ref_coord = np.genfromtxt(input.reference_coord, delimiter=',')
#        mask_cmd = "mask_data {{input.ifgram}} {wd} {{output.int_mask}} {{input.mask}} 1".format(wd=wd)
#        shell(mask_cmd)
#        filt_cmd = "fspf {{output.int_mask}} {{output.int_filt}} {wd} 0 {{params.aps_window}} 3 {{input.mli_par}}".format(wd=wd)
#        shell(filt_cmd)
#        mcf_cmd = "mcf {{output.int_filt}} - - {{output.aps_unw}} {wd} - - - - - 1 1 - {ridx} {azidx} 1 ".format(wd=wd, ridx=ref_coord[0], azidx=ref_coord[1])
#        shell(mcf_cmd)
#        #interpolate
#        interp_cmd = "interp_ad {{output.aps_unw}} {{output.aps}} {wd} 250 10 250 3 2".format(wd=wd)
#        shell(interp_cmd)
#        #Unwrap
#        mcf_cmd = "mcf {{input.ifgram}} {{input.cc}} {{input.cc_mask}} {{output.int_unw}} {wd} - - - - - 1 1 - {ridx} {azidx} 1 ".format(wd=wd, ridx=ref_coord[0], azidx=ref_coord[1])
#        shell(mcf_cmd)
#        #Mask out moving areas
#        mask_cmd = "mask_data {{output.int_unw}} {wd} {{output.int_mask}} {{input.mask}} 0".format(wd=wd)
#        shell(mask_cmd)
#        #Filter
#        filt_cmd = "fspf {{output.int_unw}} {{output.aps_masked}} {wd} 2 {{params.aps_window}} 3 {{input.mli_par}}".format(wd=wd)
#        shell(filt_cmd)
#        #Close holes
#        interp_cmd = "interp_ad {{output.aps_masked}} {{output.aps}} {wd} 250 10 250 3 2".format(wd=wd)
#        shell(interp_cmd)


#refine the aps using the stacked interferograms
rule aps_refinement:
    input:
        stack.previous_stack,
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    output:
        aps = 'diff/{mastername}_{slavename}.aps_ref',


rule cleanup_diff:
    run:
        shell('rm diff/*.int*')
        shell('rm diff/*.aps')
        shell('rm diff/*.diff*')
        shell('rm diff/*cc_mask*')
        shell('rm diff/*unw_mask*')
        shell('rm diff/*.unw')
rule cleanup_geo:
    run:
        shell('rm geo/Dom*')
        shell('rm geo/bisgletscher*')


#Remove the aps from the interferogram
rule diff_ifgram:
    output:
        diff_int = 'diff/{mastername}_{slavename}.diff',
    input:
        diff_par = 'diff/{mastername}_{slavename}.diff_par',
        aps = 'diff/{mastername}_{slavename}.aps',
        int = 'diff/{mastername}_{slavename}.unw',
        int_par = 'int/{mastername}_{slavename}.int_par',
        mli1_par = "mli/{mastername}.mli.par",
        mli2_par=  "mli/{slavename}.mli.par",
        mask = "diff/{mastername}_{slavename}.unw_mask.bmp",
    wildcard_constraints:
        mastername=slc_regex,
        slavename=slc_regex,
    run:
        #Subtract aps
        sub_cmd = "sub_phase {input.int} {input.aps} {input.diff_par} {output.diff_int} 0 0"
        shell(sub_cmd)
#        fit_cmd = "quad_fit {input.aps} {output.diff_par} - - {input.mask} - 0 "
#        shell(fit_cmd)
#        #Remove phase
#        sub_cmd = "quad_sub {input.int} {output.diff_par} {output.diff_int} 0"
#        shell(sub_cmd)
        from pyrat.diff.core import Interferogram as intgram
        import pyrat.fileutils.gpri_files as gpf
        import numpy as np
        ifgram = intgram(input.diff_par, output.diff_int, master_par = input.mli1_par, slave_par = input.mli2_par, dtype=gpf.type_mapping['FCOMPLEX'])
#        aps = gpf.gammaDataset(output.diff_par, input.aps, dtype=gpf.type_mapping['FLOAT'])
#        ifgram = np.exp(-1j * np.array(aps)) * ifgram
        ifgram.tofile(input.diff_par, output.diff_int)


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
        avg_ifgram  = 'stack/{start_dt}_stack_{nifgrams}_{chan}.diff',
        avg_ifgram_par  = 'stack/{start_dt}_stack_{nifgrams}_{chan}.diff_par',
    input:
        ifgrams = stack.all_diff,
        lut = 'geo/' + config['geocoding']['table_name']  + '.gpri_to_dem',
        dem_par = 'geo/' + config['geocoding']['table_name']  + '.dem_seg.par',
    wildcard_constraints:
        start_dt = dt_regex,
#        stop_dt = dt_regex,
    params:
        ridx = 1173,
        azidx = 112,
        ws = 8
    script:
        'scripts/stacking.py'


rule variogram:
    output:
        avg_ifgram  = 'stack/{start_dt}_stack_{nifgrams}_{chan}.variogram',
    input:
        ifgrams = stack.all_diff,
        lut = 'geo/' + config['geocoding']['table_name']  + '.gpri_to_dem',
        dem_par = 'geo/' + config['geocoding']['table_name']  + '.dem_seg.par',
    wildcard_constraints:
        start_dt = dt_regex,
#        stop_dt = dt_regex,
    params:
        ridx = 1173,
        azidx = 112,
        ws = 8
    script:
        'scripts/variogram.py'

#Here we group the rules for point target analysis


#Apply the deformation model
rule def_mod:
    output:
        aps = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pres',#residual
        dh = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pdh',#heigh correction
        defo ='ipt/{start_dt}_stack_{nifgrams}_{chan}.pdef',#deformation model
        unw ='ipt/{start_dt}_stack_{nifgrams}_{chan}.punw',#unwrapped differential phase
        sigma ='ipt/{start_dt}_stack_{nifgrams}_{chan}.psigma',#unwrapped differential phase
        pmask ='ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask_mod',#mask of accepted points
    input:
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',#list of point values
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',#list of point values
        pSLC_par ='ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC_par',#list of point values
        pitab = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pitab',
#        pbase = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pbase',
        pint =  'ipt/{start_dt}_stack_{nifgrams}_{chan}.pint',
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '.json',
        mli_par = stack.first_valid_mli_par,
        slc_par_names = (lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards))

    run:
        import pyrat.ipt.core as ipt
        import pyrat.fileutils.gpri_files as gpf
        mli_par = gpf.par_to_dict(input.mli_par)
        #Get location closest to reference
        with open(input.reference_coord) as inputfile:
            ref_coord =  get_reference_coord(json.load(inputfile))
        #load plist
        plist = ipt.Plist(input.plist, input.slc_par_names[0])
        #Convert to slc geometryd
        ref_coord = ref_coord[0] * mli_par.range_looks, ref_coord[1] * mli_par.azimuth_looks
        ref_idx = plist.closest_index(ref_coord)
        cmd = "def_mod_pt {{input.plist}} {{input.pmask}} {{input.pSLC_par}} - {{input.pitab}} - - {{input.pint}} 1 {ref_idx} {{output.aps}} {{output.dh}} {{output.defo}} {{output.unw}} {{output.sigma}} {{output.pmask}} 0 -1e3 1e3 - 5 - - - - -".format(ref_idx=ref_idx)
        shell(cmd)


rule unwrap_plist:
    output:
        aps ='ipt/{start_dt}_stack_{nifgrams}_{chan}.paps',#unwrapped differential phase
    input:
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',#list of point values
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',#list of point values
        pSLC_par ='ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC_par',#list of point values
        pitab = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pitab',
        pint =  'ipt/{start_dt}_stack_{nifgrams}_{chan}.pint',
        reference_coord = 'geo/' + config['geocoding']['table_name'] + '.json',
        mli_par = stack.first_valid_mli_par,
        slc_par_names = (lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)),
        glacier_mask = 'geo/bisgletscher_mask.bmp',
        pcc = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pccs'
    run:
        import pyrat.ipt.core as ipt
        import pyrat.fileutils.gpri_files as gpf
        mli_par = gpf.par_to_dict(input.mli_par)
        #Get location closest to reference
        with open(input.reference_coord) as inputfile:
            ref_coord =  get_reference_coord(json.load(inputfile))
        #load plist
        plist = ipt.Plist(input.plist, input.slc_par_names[0])
        #Convert to slc geometryd
        ref_coord = ref_coord[0] * mli_par.range_looks, ref_coord[1] * mli_par.azimuth_looks
        ref_idx = plist.closest_index(ref_coord)
        mcf_cmd = "mcf_pt  {{input.plist}} {{input.pmask}} {{input.pint}} - {{input.pcc}} {{input.glacier_mask}} {{output.aps}} - - {ref_idx} - {rlks} {azlks}".format(ref_idx=ref_idx, rlks=mli_par.range_looks, azlks=mli_par.azimuth_looks)
        shell(mcf_cmd)
#Determine point scatterer candidate list
rule plist:
    output:
        plist =  'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist',
        MSR =  'ipt/{start_dt}_stack_{nifgrams}_{chan}.MSR',
    input:
        slc_tab =  'ipt/{start_dt}_stack_{nifgrams}_{chan}.slc_tab',
        slc_par =  'slc_desq/{start_dt}_{chan}.slc_dec.par'
    run:
        stat_cmd  = "pwr_stat {input.slc_tab} {input.slc_par} {output.MSR} {output.plist} {config[ptarg][msr_min]} {config[ptarg][pwr_min]} - - - - 0 0"
        shell(stat_cmd)

#Mask glacier:
rule plist_mask:
    output:
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',
    input:
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist',
        mask = 'geo/bisgletscher_mask.bmp'
    run:
        cmd = "msk_pt {input.plist} - {input.mask} {output.plist} {output.pmask} {config[interferogram][rlks]} {config[interferogram][azlks]}"
        shell(cmd)


#Point list for the glacier
rule plist_glacier:
    output:
        plist_temp = temp('ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_temp'),
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_glacier',
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask_glacier',
    input:
         mask = 'geo/bisgletscher_mask_inv.bmp',
         slc_par =  'slc_desq/{start_dt}_{chan}.slc_dec.par'
    run:
        import pyrat.fileutils.gpri_files as gpf
        slc_par = gpf.par_to_dict(input.slc_par)
        cmd = "mkgrid {{output.plist_temp}} {nr} {naz} 10 10 - -".format(nr=slc_par['range_samples'], naz=slc_par['azimuth_lines'])
        shell(cmd)
        mask_cmd = "msk_pt {output.plist_temp} - {input.mask}  {output.plist} {output.pmask} {config[interferogram][rlks]} {config[interferogram][azlks]}"
        shell(mask_cmd)

#Create slc_tab:
rule slc_tab:
    output:
        slc_tab = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.slc_tab'
    input:
        slc_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec', wildcards),
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)
    run:
        with open(output.slc_tab, 'w+') as of:
            for slc, slc_par in zip(input.slc_names, input.slc_par_names):
                of.write(slc + ' '+  slc_par + '\n')

#Create pSLC_par:
rule pSLC_par:
    output:
        pSLC_par = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC_par',
        pSLC = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC'
    input:
        slc_pars = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards),#in theory we could take the slcs from the tab, but by creating it we make sure that all of them exist
        slc_tab = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.slc_tab',
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',
        slc_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec', wildcards),
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)
    run:
        cmd = 'SLC2pt {input.slc_tab} {input.plist} {input.pmask} {output.pSLC_par} {output.pSLC} -'
        shell(cmd)


#Create itab
rule pitab:
    output:
        pitab = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pitab'
    input:
        slc_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec', wildcards),
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)
    run:
        import itertools
        with open(output.pitab, 'w+') as of:
            counter = 0
            for image_counter, idx_master in enumerate(range(1, int(wildcards.nifgrams))):
                for idx_slave in range(idx_master + 1, idx_master + 2, 1):
                    print(idx_master, idx_slave)
                    if idx_slave < int(wildcards.nifgrams):
                        counter += 1
                        of.write("{idx_master} {idx_slave} {counter} 1 \n".format(idx_master=idx_master, idx_slave=idx_slave, counter=counter))




#Interferograms for point data stack
rule pint:
    output:
        pint = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pint'
    input:
        slc_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec', wildcards),
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)   ,
        pitab = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pitab',
        pSLC = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC',
        pSLC_par = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC_par',
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',
    run:
        cmd = "intf_pt {input.plist} {input.pmask} {input.pitab} - {input.pSLC} {output.pint} 0 {input.pSLC_par}"
        shell(cmd)


#Temporal Coherence for point data stack
rule pcct:
    output:
        pcc = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pcct'
    input:
        slc_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec', wildcards),
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)   ,
        pint = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pint',
        pSLC_par = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC_par',
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',
    run:
        cmd = "cct_pt {input.plist} {input.pmask} {input.slc_par_names[0]} {input.pint} {output.pcc} 0 - - -"
        shell(cmd)

#Spatial Coherence for point data stack
rule pccs:
    output:
        pcc = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pccs'
    input:
        slc_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec', wildcards),
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)   ,
        pint = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pint',
        pSLC_par = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pSLC_par',
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',
    run:
        cmd = "ccs_pt {input.plist} {input.pmask} {input.slc_par_names[0]} {input.pint} {output.pcc} - 0 - - -"
        shell(cmd)


#Reduces the point density for a point stack
rule density_reduction:
    output:
        reduced_mask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask_density_reduction',
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_reduced',
#        pdata = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.{pdata}_reduced',
    input:
        pmask = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pmask',
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards),
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist_masked',
        pquality = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.pcct',
#        pdata = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.{pdata}',
    run:
        cmd = 'pt_density_reduction {input.plist} {input.pmask} {input.slc_par_names[0]} {input.pquality} 20 {output.reduced_mask}'
        shell(cmd)
        mask_cmd = "msk_pt {input.plist} {output.reduced_mask} - {output.plist} - - -"
        shell(mask_cmd)

#Height for points in the stack
rule hgt_pt:
    input:
        plist = 'ipt/{start_dt}_stack_{nifgrams}_{chan}.plist{type}',
        hgt = 'geo/Dom.dem_seg_fgc',
        mli_par = stack.first_valid_mli_par,
        slc_par_names = lambda wildcards: stack.all_single('slc_desq/{date}_{chan}.slc_dec.par', wildcards)
    output:
        phgt =  'ipt/{start_dt}_stack_{nifgrams}_{chan}.phgt{type,(.*)|(_\w+)}',
    run:
        cmd = "data2pt {input.hgt} {input.mli_par} {input.plist} {input.slc_par_names[0]} {output.phgt} 1 2"
        shell(cmd)