#!/bin/bash
# Downloads a representative sample of all available ocean observation types
#  (SST / SSS / altimetry / insitu T&S), from publically available links
#
# Modify the start/end dates, and uncomment the desired sources below.

date_start=20180415
date_end=20180415
output_path=./data


# start the loop
date=$date_start
while [[ $date -le $date_end ]]; do
    echo $date

    ##------------------------------------------------------------
    ## SSS
    ## Note: you're best off using just the jpl_smap source
    ##------------------------------------------------------------
    bash source.jpl_smap.sh $date $output_path
    #bash source.nesdis_sss.sh smos $date $output_path
    #bash source.nesdis_sss.sh smap $date $output_path


    ##------------------------------------------------------------
    ## ADT
    ##------------------------------------------------------------
    bash source.nesdis_adt_rads.sh $date $output_path
    

    ##------------------------------------------------------------
    ## insitu
    ##------------------------------------------------------------
    bash source.fnmoc.sh prof $date $output_path
    bash source.fnmoc.sh sfc  $date $output_path
    bash source.fnmoc.sh trak $date $output_path
    

    ##------------------------------------------------------------
    ## SST
    ## Note: best off using just VIIRS L3U, unless you have
    ##  a specific need for any of the others.
    ##  (L2 data is huge... especially for VIIRS)
    ##------------------------------------------------------------
    #bash source.nesdis_sst_viirs.sh l2p $date $output_path
    bash source.nesdis_sst_viirs.sh l3u $date $output_path
    #bash source.nesdis_sst_avhrr.sh l2p $date $output_path
    #bash source.nesdis_sst_avhrr.sh l3u $date $output_path
    
date=$(date -d "$date + 1 day" +%Y%m%d)
    
done
