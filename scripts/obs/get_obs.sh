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
    ## Note: best off using just the RSS 70km L2 for now
    ##------------------------------------------------------------
    bash source.jpl_smap.sh rss_70km $date $output_path
    #bash source.jpl_smap.sh rss_40km $date $output_path
    #bash source.jpl_smap.sh jpl $date $output_path


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
    ## ice concentration
    ##------------------------------------------------------------
    bash source.nsidc_icec.sh $date $output_path

    ##------------------------------------------------------------
    ## ice freeboard/thickness
    ##------------------------------------------------------------
    bash source.esa_ice_cryosat.sh GDR $date $output_path
    #bash source.esa_ice_cryosat.sh LRM $date $output_path
    #bash source.esa_ice_cryosat.sh SAR $date $output_path
    #bash source.esa_ice_cryosat.sh SIN $date $output_path

    ##------------------------------------------------------------
    ## SST
    ## Note: best off using just VIIRS/AVHRR L3U, and the microwave
    ## obs, unless you have a specific need for any of the others.
    ##  (L2 data is huge... especially for VIIRS)
    ##------------------------------------------------------------
    bash source.ghrsst_sst.sh amsr2 l3u $date $output_path
    bash source.ghrsst_sst.sh gmi l3u $date $output_path
    #bash source.ghrsst_sst.sh goes16 l3u $date $output_path
    bash source.ghrsst_sst.sh windsat l3u $date $output_path
    #bash source.nesdis_sst_viirs.sh l2p $date $output_path
    bash source.nesdis_sst_viirs.sh l3u $date $output_path
    #bash source.nesdis_sst_avhrr.sh l2p $date $output_path
    bash source.nesdis_sst_avhrr.sh l3u $date $output_path

    ##------------------------------------------------------------
    ## Ocean Color
    ## L2 data from VIIRS/JPSS1, VIIRS/SNPP, MODIS-Aqua
    ##------------------------------------------------------------
    #bash source.coastwatch_oc_viirs_modis.sh snpp $date $output_path
    bash source.coastwatch_oc_viirs_modis.sh jpss1 $date $output_path
    bash source.coastwatch_oc_viirs_modis.sh aqua $date $output_path 

date=$(date -d "$date + 1 day" +%Y%m%d)

done
