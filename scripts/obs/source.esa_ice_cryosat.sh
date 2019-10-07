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
if [[ $radar == "GDR" ]]; then
   source="ftp://science-pds.cryosat.esa.int/SIR_${radar}/${yr}/${mm}/"
else
   source="ftp://science-pds.cryosat.esa.int/SIR_${radar}_L2/${yr}/${mm}/"
fi

pwd=$(pwd)
d=$out_dir/$date
mkdir -p $d
cd $d

files=$(curl $source -l)
namelist1=*__$date*.DBL
namelist2=*__$date*.nc

for f in $files; do
    # Only keep .dbl or .nc files
    if [[ $f == $namelist1 ]] || [[ $f == $namelist2 ]]; then
      wget $source/$f 
    fi 
done
cd $pwd
