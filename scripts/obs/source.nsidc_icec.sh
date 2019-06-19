#!/bin/bash
# download ice concentration observations from NSIDC.

set -e
set -u

usage="usage: $0 yyyymmdd output_path"
if [[ $# != 2 ]]; then
    echo $usage
    exit 1
fi

date=$1
yr=${date:0:4}
mm=${date:4:2}
dy=${date:6:2}

out_dir="$2/icec.nsidc"
sourcenh="https://www.ncei.noaa.gov/data/sea-ice-concentration/access/${yr}/seaice_conc_daily_nh_f17_${yr}${mm}${dy}_v02r00.nc"
sourcesh="https://www.ncei.noaa.gov/data/sea-ice-concentration/access/${yr}/seaice_conc_daily_sh_f17_${yr}${mm}${dy}_v02r00.nc"

pwd=$(pwd)
d=$out_dir/$date
mkdir -p $d
cd $d
wget $sourcenh
wget $sourcesh
cd $pwd
