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


source_base="https://podaac-tools.jpl.nasa.gov/drive/files/allData/smap/L2"

if [[ $type == "jpl" ]]; then
    source=$source_base"/JPL/V4.2"
    file_sfx='h5'
elif [[ $type == "rss" ]]; then
    source=$source_base"/RSS/V4/SCI"
    file_sfx='nc'
else
    echo $usage
    exit 1
fi

out_dir="$output_path/sss.smap.$type"
pwd=$(pwd)
source_dir=$source/${yr}/${dy}/

d=$out_dir/$date
mkdir -p $d
cd $d

wget -nd -r -A $file_sfx --user=myusername --password=mypassword $source_dir

cd $pwd
