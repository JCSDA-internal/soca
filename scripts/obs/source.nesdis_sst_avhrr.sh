#!/bin/bash
# Downloads NESDIS SST ACSPO retrievals for AVHRR L2/L3
set -e
set -u

usage="usage: $0 [l2p|l3u] yyyymmdd output_path"
if [[ $# != 3 ]]; then
    echo $usage
    exit 1
fi

date=$2
yr=${date:0:4}
dy=$(date -d "$date" "+%j")


lvl=$1
valid=0
for l in l2p l3u; do
    [[ $l == $lvl ]] && valid=1
done
if [[ $valid != 1 ]]; then
    echo $usage
    exit 1
fi

out_dir="$3/sst.avhrr_${lvl}.nesdis"
source="ftp://ftp.star.nesdis.noaa.gov/pub/socd2/coastwatch/sst/ran/avhrr_gac/"
sats="metopa noaa16 noaa17 noaa18 noaa19"

pwd=$(pwd)

for sat in $sats; do
    source_dir=$source/$sat/$lvl/$yr/$dy/
    echo $source_dir
    files=$(curl -lf $source_dir || echo "" )
    for f in $files; do
	d=$out_dir/$date
	mkdir -p $d
	cd $d
	wget $source_dir/$f
	cd $pwd
    done
    
done

