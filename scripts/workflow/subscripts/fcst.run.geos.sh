#!/bin/bash
set -e
cat <<EOF



#================================================================================
#================================================================================
# fcst.run.sh
#  Runs a MOM6/SIS2 forecast
#================================================================================

EOF

# required environment variables:
envar=()
envar+=("FCST_LEN")          # Length of the forecast (hours)
envar+=("FCST_RESTART")      # =1 if a restart file is used, otherwise T/S IC file is used
envar+=("FCST_START_TIME")   # start of forecast (in any appropriate "date" command format)
envar+=("MOM_CONFIG")        # path to input model configuration files
envar+=("MOM_DATA")          # path to input model static data
envar+=("MOM_EXE")           # path to MOM6/SIS2 executable
envar+=("MOM_IC")            # path to T/S IC. Only used if FCST_RESTART==0
envar+=("RESTART_DIR_IN")    # path to restart files from previous cycle (if FCST_RESTART==1)
envar+=("WORK_DIR")          # temporary working directory for this script

# make sure required env vars exist
set +u
for v in ${envar[@]}; do
    if [[ -z "${!v}" ]]; then
    echo "ERROR: env var $v is not set."; exit 1
    fi
    echo " $v = ${!v}"
done
set -u
echo -e "\n================================================================================\n"

# ================================================================================

FCST_RST_OFST=$FCST_LEN


#--------------------------------------------------------------------------------

# create working directory
rm -rf $WORK_DIR
mkdir $WORK_DIR
cd $WORK_DIR
mkdir OUTPUT
mkdir RESTART

geos_scripts=/home/gvernier/nobackup/geos-santha/ #TODO(Guillaume) move upstream

# prepare resource files for geos
cp -r /home/gvernier/nobackup/geos-santha/geos-scratch/* $WORK_DIR

# generate cap_restart
cap_restart=$(date --utc "+%Y%m%d %H%M%S" -d "$FCST_START_TIME")
cat <<EOF > cap_restart 
$(date --utc "+%Y%m%d %H%M%S" -d "$FCST_START_TIME")
EOF

# link mom6 files
ln -s $MOM_CONFIG/* .

# link model static input files
mkdir INPUT
cd INPUT
ln -s $MOM_DATA/* .

# get ogcm ICs
if [[ $FCST_RESTART == 1 ]]; then
    ln -s $RESTART_DIR_IN/MOM*.nc .
else
    ln -s $MOM_IC ic.nc
fi
cd ../

# get agcm ICs
if [[ $FCST_RESTART == 1 ]]; then
    cp $RESTART_DIR_IN/*_rst .
else
    #TODO(Guillaume) make location of initial agcm restarts user defined
    date_str=$(date --utc "+%Y%m%d_%H%M%S" -d "$FCST_START_TIME")
    /bin/cp /home/gvernier/nobackup/geos-santha/soca-cpld/util/rst_atm_${date_str}.tar $WORK_DIR
    tar -xvf rst_atm_${date_str}.tar
fi

# create the time dependent mom6 namelist file
. input.nml.sh > input.nml

# run geos
$geos_scripts/gcm_run.j

# prepare aogcm restart for next cycle
mv ./scratch/RESTART/* ./RESTART/
mv ./scratch/*_checkpoint ./RESTART/

# rename agcm restarts
cd ./RESTART/
checkpoints=`/bin/ls -1 *_checkpoint`
for checkpoint in $checkpoints
do
  rename _checkpoint _rst ${checkpoint}
done
