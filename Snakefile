import pyrat
from snakemake.remote.SFTP import RemoteProvider, RemoteObject
import glob
provider = RemoteProvider(username="baffelli", private_key="/home/baffelli/.ssh/id_rsa")
data_source = RemoteObject(provider=provider)

server = "ifu-eo-srv-1.ethz.ch/local/unique_data/2015_GPRI_Dom/"

#subworkflow new_data:
#    workdir: './processed'
#    snakefile:  pyrat.rules['slc_to_calibrated_c']
#    configfile: './bisgletscher.json'
#
date = "20150710_1*"
print(glob.glob(provider.remote(server + 'raw/*')))


rule all:
    input:
        expand('slc_chan/20150710_174104_AAAl.slc'),
        'slc_chan/20150710_174334_AAAl.slc'

rule untar_and_copy:
        input:
            tar = provider.remote(server + "raw/{datetime}.tar", static=True),
        output:
            slc = 'slc_chan/{datetime}_{chan}.slc',
            slc_par = 'slc_chan/{datetime}_{chan}.slc.par'
        run:
            shell("tar -xvf {input.tar} --directory ./slc_chan --strip=3 data/Dom/slc/{wildcards.datetime}_{wildcards.chan}.slc")
            shell("tar -xvf {input.tar} --directory ./slc_chan --strip=3 data/Dom/slc/{wildcards.datetime}_{wildcards.chan}.slc.par ")
