import pyrat
from snakemake.remote.SFTP import RemoteProvider

SFTP = RemoteProvider(username="baffelli", private_key="/home/baffelli/.ssh/id_rsa")


#subworkflow new_data:
#    workdir: './processed'
#    snakefile:  pyrat.rules['slc_to_calibrated_c']
#    configfile: './bisgletscher.json'
#
rule all:
    input:
        'slc_chan/20150710_174104_AAAl.slc'

rule untar_and_copy:
        input:
            tar = SFTP.remote("ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/raw/{datetime}.tar"),
        output:
            slc = 'slc_chan/{datetime}_{chan}.slc',
            slc_par ='slc_chan/{datetime}_{chan}.slc.par'
        run:
            shell("tar -xvf {input.tar} --strip=3 {wildcards.datetime}_{wildcards.chan}.slc -C {output.slc}")
            shell("tar -xvf {input.tar} --strip=3 {output.slc_par}")