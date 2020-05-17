#!/bin/bash
set -e
cat <<EOF

#================================================================================
#================================================================================
# da.run.sh
#  Run the SOCA data assimilation (currently only 3dvar non-fgat)
#================================================================================

EOF

# Required environment variables:
envar+=()
envar+=("ANA_TIME")        # date and time of analysis (in any valid "date" command format)
envar+=("CYCLE_RST_DIR")   # path to output restart directory for running next forecast
envar+=("DA_INIT_DIR")     # path to soca bkgerr/geometry files
envar+=("FCST_LEN")        # length of forecast (hours)
envar+=("FCST_START_TIME") # date and time of start of forecast
envar+=("MOM_CONFIG")      # path to input model configuration files
envar+=("MOM_DATA")        # path to input model static data
envar+=("OBS_IODA")        # path to observations already processed into ioda format
envar+=("OBS_ADT")         # directory name for adt
envar+=("OBS_INSITU")      # directory name for insitu
envar+=("OBS_SST")         # directory name for sst
envar+=("RESTART_DIR")     # path to input restart files for da background
envar+=("SOCA_BIN_DIR")    # path to soca executables
envar+=("SOCA_CONFIG")     # path to input soca configuration files
envar+=("SOCA_DATA")       # path to input soca static data
envar+=("WORK_DIR")        # temporary working directory for this script

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


ulimit -s unlimited
export OMP_NUM_THREADS=1


# create working directory
rm -rf $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR
mkdir RESTART
mkdir OUTPUT
mkdir obs
mkdir obs_out
mkdir Data

# link mom6 static files
ln -s $MOM_CONFIG/* .
mkdir INPUT
cd INPUT
ln -s $MOM_DATA/* .
cd ..


# create time dependent mom namelist file
FCST_RESTART=1
FCST_RST_OFST=$FCST_LEN
. input.nml.sh > input.nml

# link soca config files
ln -s $SOCA_CONFIG/* .
. 3dvar.yml.sh > 3dvar.yml
. checkpoint.yml.sh > checkpoint.yml
ln -s $SOCA_DATA/* .

# link bump / gridgen files
ln -s $DA_INIT_DIR/* .

# link the input restart files
cd INPUT
ln -s $RESTART_DIR/* .
cd ..

# link the observations
# TODO make obs list configurable
cd obs
ymd=$(date -u -d "$ANA_TIME" +%Y%m%d)
ln -s $OBS_IODA/$OBS_ADT/${ymd:0:4}/$ymd.nc adt.nc
ln -s $OBS_IODA/$OBS_INSITU/${ymd:0:4}/$ymd.nc insitu.nc
ln -s $OBS_IODA/$OBS_SST/${ymd:0:4}/$ymd.nc sst.nc
cd ..

# run the 3dvar
echo "Starting 3dvar..."
time mpirun -np 360 $SOCA_BIN_DIR/soca_3dvar.x 3dvar.yml

# do the checkpointing
echo "Starting checkpointing..."
time mpirun -np 360 $SOCA_BIN_DIR/soca_checkpoint_model.x checkpoint.yml

# move the restart files to a non-scratch space
echo "Moving restart files..."
mkdir -p $CYCLE_RST_DIR
cp $WORK_DIR/INPUT/*.res* $CYCLE_RST_DIR/
cp $WORK_DIR/INPUT/*_rst $CYCLE_RST_DIR/
cp $WORK_DIR/RESTART/* $CYCLE_RST_DIR/
