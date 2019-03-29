#!/bin/bash
# download sea surface salinity observations from SMAP
# TODO: in order to handle all available dates,
#  scripts needs to be changed to correctly select the correct source directory

set -e
set -u

usage="usage: $0 [jpl|rss_40km|rss_70km] yyyymmdd output_path"
if [[ $# != 3 ]]; then
    echo $usage
    exit 1
fi

type=$1
date=$2
output_path=$3


yr=${date:0:4}
dy=$(date -d "$date" "+%j")


if [[ $type == "jpl" ]]; then
    source="ftp://podaac-ftp.jpl.nasa.gov/allData/smap/L2/JPL/V4.2"
elif [[ $type == "rss_40km" ]]; then
    source="ftp://podaac-ftp.jpl.nasa.gov/allData/smap/L2/RSS/V3/SCI/40KM"
elif [[ $type == "rss_70km" ]]; then
    source="ftp://podaac-ftp.jpl.nasa.gov/allData/smap/L2/RSS/V3/SCI/70KM"
else
    echo $usage
    exit 1
fi
    
out_dir="$output_path/sss.smap.$type"
pwd=$(pwd)
source_dir=$source/${yr}/${dy}/


files=$(curl $source_dir -l)
for f in $files; do
    # ignore the .md5 files
    fe=${f##*.}
    [[  $fe == "md5" ]] && continue
    
    d=$out_dir/$date
    mkdir -p $d
    cd $d
    wget $source_dir/$f
    cd $pwd
done
