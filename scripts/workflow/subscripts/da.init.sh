#!/bin/bash
set -e
cat <<EOF

#================================================================================
#================================================================================
# da.init.sh
#  Initialize the static B if it has not already been initialized.
#  Also initialize the geometry file if not already created.
#  If these have already been created, this script will silently exit.
#================================================================================

EOF

# Required environment variables:
envar=()
envar+=("DA_INIT_DIR")     # path for output bkgerr/geometry files
envar+=("FCST_LEN")        # length of forecast (hours)
envar+=("FCST_START_TIME") # datetime of start of forecast (in any valid "date" command format)
envar+=("MOM_CONFIG")      # path to input model configuration files
envar+=("MOM_DATA")        # path to input model static data
envar+=("RESTART_DIR")     # path to input restart files for DA background
envar+=("SOCA_BIN_DIR")    # path to soca executables directory
envar+=("SOCA_CONFIG")     # path to input soca configuration files
envar+=("SOCA_DATA")       # path to input soca static data
envar+=("WORK_DIR")        # temporary working directory for this script
envar+=("MPIRUN")          # exec to run mpi
envar+=("JOB_NPES")        # exec to run mpi

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

# if static b already exists, stop here
if [[ -e "$DA_INIT_DIR/valid" ]]; then
    echo "Static B has already been initialized."
    exit 0
fi

# create working directory
rm -rf $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR
mkdir bump
mkdir OUTPUT

# link model files
ln -s $MOM_CONFIG/* .

# link model static input files
mkdir INPUT
cd INPUT
ln -s $MOM_DATA/* .
cd ..

# create time dependent momt namelist file
FCST_RESTART=1
FCST_RST_OFST=$FCST_LEN
. input.nml.sh > input-mom6.nml

# link soca config files
ln -s $SOCA_CONFIG/* .
ln -s $SOCA_DATA/* .

# link the input restart files
cd INPUT
ln -s $RESTART_DIR/* .
cd ..

# run the grid generation
ulimit -s unlimited
export OMP_NUM_THREADS=1
echo "generating grid..."
time $MPIRUN -np $JOB_NPES $SOCA_BIN_DIR/soca_gridgen.x gridgen.yml

# run static B init
echo "generating static b..."
time $MPIRUN -np $JOB_NPES $SOCA_BIN_DIR/soca_staticbinit.x staticbinit.yml

# move output files
rm -rf $DA_INIT_DIR
mkdir -p $DA_INIT_DIR
mv bump $DA_INIT_DIR/
mv soca_gridspec.nc $DA_INIT_DIR/
touch $DA_INIT_DIR/valid
