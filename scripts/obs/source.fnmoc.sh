#!/bin/bash
# Downloads insitu ocean observations from FNMOC (GODAE)

set -e
set -u

usage="usage: $0 [prof|sfc|trak] yyyymmdd output_path"
if [[ $# != 3 ]]; then
    echo $usage
    exit 1
fi

date=$2
yr=${date:0:4}
dy=$(date -d "$date" "+%j")

type=$1

out_dir="$3/insitu.fnmoc"
source="https://www.usgodae.org/pub/outgoing/fnmoc/data/ocn"


if [[ $type == "prof" ]]; then
    fn1=profile
    fn2=profile
elif [[ $type == "sfc" ]]; then
    fn1=sfcobs
    fn2=ship
elif [[ $type == "trak" ]]; then
    fn1=trak
    fn2=trak
else
    echo $usage
    exit 1
fi

d=$out_dir/$date
pwd=$(pwd)
mkdir -p $d
cd $d
wget $source/$fn1/$yr/${date}00.$fn2.Z
gunzip *.Z
cd $pwd
