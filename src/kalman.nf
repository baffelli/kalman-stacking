#! ../nextflow


all_raws = Channel.fromPath( "/mnt/unique_data/2015_GPRI_Dom/raw/20150803*.tar" )


process extract_slc {

    output:
    file 'slc' into slc_chan
    file 'slc_par' into slc_chan_par

    input:
    file tar from all_raws

    """
        tar -xvf $tar --directory ./slc_chan --strip=3 $slc_chan $slc_chan_par
    """

}


