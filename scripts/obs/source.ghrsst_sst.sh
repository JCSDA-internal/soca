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
source_base="ftp://podaac-ftp.jpl.nasa.gov/allData/ghrsst/data/GDS2"
ignore="^$"
if [[ $sat == "amsr2" ]]; then
    source="$source_base/$lvlU/AMSR2/REMSS/v8a/"
    ignore="rt-"
elif [[ $sat == "gmi" ]]; then
    source="$source_base/$lvlU/GMI/REMSS/v8.2a/"
elif [[ $sat == "goes16" ]]; then
    source="$source_base/$lvlU/GOES16/OSPO/v2.5/"
elif [[ $sat == "windsat" ]]; then
    source="$source_base/$lvlU/WindSat/REMSS/v7.0.1a/"
    ignore="rt-"
else
    echo $usage
    exit 1
fi 

out_dir="$4/sst.${sat}_${lvl}.ghrsst"

pwd=$(pwd)

source_dir=$source/$yr/$dy/
echo $source_dir
files=$(curl -lf $source_dir || echo "" )
for f in $files; do
    # ignore the .md5 files
    fe=${f##*.}
    [[ $fe == "md5" ]] && continue

    # some other files need to be ignored, depending on the source
    [[ $f =~ $ignore ]] && continue
    
    d=$out_dir/$date
    
    # skip if file already exists
    [[ -e $d/$f ]] && continue

    mkdir -p $d
    cd $d
    wget $source_dir/$f
    cd $pwd
done

