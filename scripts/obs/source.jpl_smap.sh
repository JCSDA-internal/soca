#!/bin/bash
# download sea surface salinity observations from SMAP
# TODO: in order to handle all available dates,
#  scripts needs to be changed to correctly select the correct source directory

set -e
set -u

if [[ $# != 2 ]]; then
    echo "usage: $0 yyyymmdd output_path"
    exit 1
fi

date=$1
yr=${date:0:4}
dy=$(date -d "$date" "+%j")

out_dir="$2/sss.smap.jpl"
pwd=$(pwd)
source="ftp://podaac-ftp.jpl.nasa.gov/allData/smap/L2/JPL/V4.2"
source_dir=$source/${yr}/${dy}/


files=$(curl $source_dir -l)
for f in $files; do
    # ignore the .md5 files
    fe=${f##*.}
    [[ ! $fe == "h5" ]] && continue
    
    d=$out_dir/$date
    mkdir -p $d
    cd $d
    wget $source_dir/$f
    cd $pwd
done
