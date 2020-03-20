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
envar+=("FORC_DIR")          # path to input atmospheric forcing files
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
echo ""

# ================================================================================

FCST_RST_OFST=$FCST_LEN


#--------------------------------------------------------------------------------


ulimit -s unlimited
export OMP_NUM_THREADS=1


# create working directory
rm -rf $WORK_DIR
mkdir $WORK_DIR
cd $WORK_DIR
mkdir OUTPUT
mkdir RESTART

# link model files
ln -s $MOM_CONFIG/* .
ln -s $MOM_EXE .

if [[ "${FORC_ATM}" == "cfsr" ]]; then
    ln -sf data_table.cfsr data_table
elif [[ "${FORC_ATM}" == "merra2" ]]; then
    ln -sf data_table.merra2 data_table
fi

# link model static input files and ICs
mkdir INPUT
cd INPUT
ln -s $MOM_DATA/* .
if [[ $FCST_RESTART == 1 ]]; then
    ln -s $RESTART_DIR_IN/* .
else
    ln -s $MOM_IC ic.nc
fi
cd ../

# link the forcing
ln -s $FORC_DIR FORC

# create the time dependent mom6 namelist file
. input.nml.sh > input.nml

# run it!
mpirun ./MOM6
