#!/bin/bash
# download altimetry data from NESDIS.

set -e
set -u

if [[ $# != 2 ]]; then
    echo "usage: $0 yyyymmdd output_path"
    exit 1
fi

date=$1
yr=${date:0:4}
dy=$(date -d "$date" "+%j")

out_dir="$2/adt.nesdis"
source="ftp://ftp.star.nesdis.noaa.gov/pub/sod/lsa/rads/adt/${yr}/"

pwd=$(pwd)
files=$(curl $source -l)
for f in $files; do
    # make sure it is the right day
    [[ ! $f =~ ^rads_adt_.._$yr$dy.*$ ]] && continue

    d=$out_dir/$date
    mkdir -p $d
    cd $d
    wget $source/$f
    cd $pwd
done
