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
envar+=("FCST_LEN")
envar+=("FCST_RESTART")
envar+=("FCST_START_TIME")
envar+=("FORC_DIR")
envar+=("MOM_CONFIG")
envar+=("MOM_DATA")
envar+=("MOM_EXE")
envar+=("MOM_IC")
envar+=("WORK_DIR")
envar+=("RESTART_DIR_IN")

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
