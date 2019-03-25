#!/bin/bash
# download sea surface salinity observations from NESDIS.
#  Either SMAP or SMOS.

set -e
set -u

usage="usage: $0 yyyymmdd [smos|smap] output_path"
if [[ $# != 3 ]]; then
    echo $usage
    exit 1
fi

date=$2
yr=${date:0:4}
dy=$(date -d "$date" "+%j")

sat=$1

out_dir="$3/sss.$sat.nesdis"
if [[ $sat == smos ]]; then
    source="https://www.star.nesdis.noaa.gov/data/socd1/coastwatch/products/miras/nc/SM_D${yr}${dy}_Map_SATSSS_data_1day.nc"
elif [[ $sat == smap ]]; then
    source="https://www.star.nesdis.noaa.gov/data/socd1/coastwatch/products/smap/nc/SP_D${yr}${dy}_Map_SATSSS_data_1day.nc"
else
    echo $usage
    exit 1
fi

pwd=$(pwd)
d=$out_dir/$date
mkdir -p $d
cd $d
wget $source
cd $pwd
