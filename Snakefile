import pyrat
from snakemake.remote.SFTP import RemoteProvider, RemoteObject
import glob
provider = RemoteProvider(username="baffelli", private_key="/home/baffelli/.ssh/id_rsa")
server = "ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/raw/"
tar = provider.glob_wildcards(server + "{date,20150803}_{time,06[0-9]{4}}.tar")

configfile: './bisgletscher.json'



all_ifgrams = expand("diff/{dte}_{time}_AAAl_{dte_1}_{time_1}_AAAl.int",zip, dte=tar.date, dte_1=tar.date, time=tar.time, time_1=tar.time[1:])


include: pyrat.rules['raw_to_slc']


rule all:
    input:
       all_ifgrams

#Do not need to perform RC if the data comes from the server
ruleorder: untar_and_copy > range_compression
rule untar_and_copy:
        input:
            tar = provider.remote(server + "{datetime}.tar", static=True),
        output:
            slc = 'slc_chan/{datetime}_{chan}.slc',
            slc_par = 'slc_chan/{datetime}_{chan}.slc.par'
        run:
            shell("tar -xvf {input.tar} --directory ./slc_chan --strip=3 data/Dom/slc/{wildcards.datetime}_{wildcards.chan}.slc")
            shell("tar -xvf {input.tar} --directory ./slc_chan --strip=3 data/Dom/slc/{wildcards.datetime}_{wildcards.chan}.slc.par ")

rule ifgram:
    input:
        master = 'slc_desq/{mastername}.slc_dec',
        master_par = 'slc_desq/{mastername}.slc_dec.par',
        slave = 'slc_desq/{slavename}.slc_dec',
        slave_par = 'slc_desq/{slavename}.slc_dec.par',
    output:
        off_par = 'diff/{mastername}_{slavename}.off_par',
        ifgram = 'diff/{mastername}_{slavename}.int',
    params:
        mli1 = 'mli/{mastername}.mli',
        mli2 = 'mli/{slavename}.mli'
    run:
        shell("create_offset {input.master_par} {input.slave_par} {output.off_par} - 2 2 0")
        shell("interf_SLC {input.master} {input.slave} {input.master_par} {input.slave_par} {input.off_par} {params.mli1} {params.mli2} {output.ifgram} 2 2 - - 0 0")