#!/bin/bash
# download sea ice data from earth.esa

set -e
set -u

if [[ $# != 3 ]]; then
    echo "usage: $0 radar yyyymmdd output_path"
    exit 1
fi

radar=$1
date=$2
yr=${date:0:4}
mm=${date:4:2}
dy=${date:6:2}

out_dir="$3/ice.cryosat_${radar}.esa"
source="ftp://science-pds.cryosat.esa.int/SIR_${radar}_L2/${yr}/${mm}"

pwd=$(pwd)
d=$out_dir/$date
mkdir -p $d
cd $d

files=$(curl $source -l)
for f in $files; do
    # make sure it is the right day
    [[ ! $f =~ ^CS*$yr$mm$dy*$ ]] && continue
    wget $source/$f
done
cd $pwd
