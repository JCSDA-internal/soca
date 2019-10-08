#!/bin/bash
set -e
cat <<EOF

#================================================================================
#================================================================================
# da.init.sh
#  Initialize the static B if it has not already been initialized
#================================================================================

EOF

# Required environment variables:
envar=()
envar+=("DA_INIT_DIR")
envar+=("FCST_LEN")
envar+=("FCST_START_TIME")
envar+=("MOM_CONFIG")
envar+=("RESTART_DIR")
envar+=("SOCA_BIN_DIR")
envar+=("SOCA_CONFIG")
envar+=("WORK_DIR")

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
rm INPUT
mkdir INPUT
cd INPUT
ln -s $MOM_CONFIG/INPUT/* .
cd ..

# create time dependent momt namelist file
FCST_RESTART=1
FCST_RST_OFST=$FCST_LEN
. input.nml.sh > input.nml

# link soca config files
ln -s $SOCA_CONFIG/* .

# link the input restart files
cd INPUT
ln -s $RESTART_DIR/* .
cd ..

# run the grid generation
ulimit -s unlimited
export OMP_NUM_THREADS=1
echo "generating grid..."
time mpirun $SOCA_BIN_DIR/soca_gridgen.x gridgen.yml

# run static B init
echo "generating static b..."
time mpirun $SOCA_BIN_DIR/soca_staticbinit.x staticbinit.yml

# move output files
rm -rf $DA_INIT_DIR
mkdir -p $DA_INIT_DIR
mv bump $DA_INIT_DIR/
mv soca_gridspec.nc $DA_INIT_DIR/
touch $DA_INIT_DIR/valid
