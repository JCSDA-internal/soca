#!/bin/bash
# Downloads NESDIS SST ACSPO retrievals for VIIRS L2/L3
set -e
set -u

usage="usage: $0 [lsp|l3u] yyyymmdd output_path"
if [[ $# != 3 ]]; then
    echo $usage
    exit 1
fi

date=$2
yr=${date:0:4}
dy=$(date -d "$date" "+%j")

lvl=$1
if [[ $lvl == "l2p" ]]; then
    source="ftp://podaac-ftp.jpl.nasa.gov/allData/ghrsst/data/GDS2/L2P/VIIRS_NPP/OSPO/v2.41"
elif [[ $lvl == "l3u" ]]; then    
    source="ftp://ftp.star.nesdis.noaa.gov/pub/socd2/coastwatch/sst/ran/viirs/npp/l3u"
else
    echo $usage
    exit 1
fi

out_dir="$3/sst.viirs_${lvl}.nesdis"

source_dir=$source/${yr}/${dy}/
pwd=$(pwd)

files=$(curl $source_dir -l)
for f in $files; do
    # ignore the .md5 files
    fe=${f##*.}
    [[ ! $fe == "nc" ]] && continue
    
    d=$out_dir/$date
    mkdir -p $d
    cd $d
    wget $source_dir/$f
    cd $pwd
done
