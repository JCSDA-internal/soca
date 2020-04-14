#!/bin/bash
# Downloads NOAA Coastwatch retrievals for L2 Ocean Color
# MODIS-AQUA, VIIRS-JPSS1(NOAA-20), VIIRS-SNPP

set -e
set -u

usage="usage: $0 [jpss1|snpp|aqua] yyyymmdd output_path"
if [[ $# != 3 ]]; then
    echo $usage
    exit 1
fi

date=$2
yr=${date:0:4}
dy=$(date -d "$date" "+%j")

lvl=$1
if [[ $lvl == "snpp" ]]; then
    source="ftp://ftpcoastwatch.noaa.gov/pub/socd1/mecb/coastwatch/viirs/nrt/L2"
elif [[ $lvl == "jpss1" ]]; then    
    source="ftp://ftpcoastwatch.noaa.gov/pub/socd2/mecb/coastwatch/viirs/n20/nrt/L2"
elif [[ $lvl == "aqua" ]]; then
    source="https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/L2"
else
    echo $usage
    exit 1
fi

file_sfx='*.nc'

if [[ $lvl == "aqua" ]]; then
    out_dir="$3/oc.modis_${lvl}.nasagsfc"
else
    out_dir="$3/oc.viirs_${lvl}.coastwatch"
fi

source_dir=$source/${yr}/${dy}/
echo $source_dir
pwd=$(pwd)
d=$out_dir/$date
mkdir -p $d
cd $d

if [[ $lvl == "aqua" ]]; then
    wget -q -O - $source_dir |grep OC| wget --user=USERNAME --password=PASSWD --auth-no-challenge=on --base https://oceandata.sci.gsfc.nasa.gov/ -N --wait=0.5 --random-wait --force-html -i -
else
    wget -r -nc -np -nH -nd -A $file_sfx  $source_dir
fi
    
cd $pwd

