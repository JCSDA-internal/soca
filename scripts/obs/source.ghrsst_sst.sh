#!/bin/bash
# Downloads SST retrievals from the various GHRSST groups
set -e
set -u

usage="usage: $0 [amsr2|gmi|goes16|windsat] [l2p|l3u] yyyymmdd output_path"
if [[ $# != 4 ]]; then
    echo $usage
    exit 1
fi


date=$3
yr=${date:0:4}
dy=$(date -d "$date" "+%j")


lvl=$2
if [[ $lvl == "l2p" ]]; then
    lvlU="L2P"
elif [[ $lvl == "l3u" ]]; then
    lvlU="L3U"
else
    echo $usage
    exit 1
fi

sat=$1
source_base="https://podaac-tools.jpl.nasa.gov/drive/files/allData/ghrsst/data/GDS2"
ignore="^$"
file_sfx='*.nc'

if [[ $sat == "amsr2" ]]; then
    source="$source_base/$lvlU/AMSR2/REMSS/v8a"
elif [[ $sat == "gmi" ]]; then
    source="$source_base/$lvlU/GMI/REMSS/v8.2a"
elif [[ $sat == "goes16" ]]; then
    source="$source_base/$lvlU/GOES16/OSPO/v2.5"
elif [[ $sat == "windsat" ]]; then
    source="$source_base/$lvlU/WindSat/REMSS/v7.0.1a"
else
    echo $usage
    exit 1
fi 

out_dir="$4/sst.${sat}_${lvl}.ghrsst"

pwd=$(pwd)

source_dir=$source/$yr/$dy/
echo $source_dir

d=$out_dir/$date
mkdir -p $d
cd $d

wget -r -nc -np -nH -nd -A $file_sfx  $source_dir

cd $pwd
