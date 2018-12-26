#!/bin/sh

#source anaconda2/Python/bin/activate nco

# SOCA Initialization
#SOCA_DA_IC=false
SOCA_DA_IC=true

# Singularity image to use
JCSDA_SING="/home/gvernier/Sandboxes/soca/JCSDA-singularity-master-latest.simg"
JCSDA_SING_old="/home/gvernier/Sandboxes/soca/JCSDA-singularity-master-latest.simg-old"

# Source, config and build location of soca
SOCA_SRC="/home/gvernier/Sandboxes/soca/soca-bundle-mom6/soca"
SOCA_BUILD="/home/gvernier/Sandboxes/soca/soca-bundle-mom6/build"
SOCA_BIN="${SOCA_BUILD}/bin"
SOCA_INPUT="${SOCA_BUILD}/soca/test/testinput"
SOCA_TEST="${SOCA_BUILD}/soca/test"
SOCA_MODEL_RSC="${SOCA_SRC}/test/Data/360x210x63"
SOCA_MODEL_FORCING="${SOCA_SRC}/test/Data/360x210x63"

# MOM6 scratch/experiment template
#MOM6_MINSCRATCH="/home/gvernier/Sandboxes/soca/mom6-scratch"
MOM6_MINSCRATCH="/home/gvernier/Sandboxes/MOM6-examples/ice_ocean_SIS2/SIS2/scratch/minscratch"

echo $SOCA_DA_IC

cd ${SOCA_TEST}

# Setup background
SOCA_BKG="/home/gvernier/Sandboxes/soca/data/mom6-rst/MOM.res.mom6-sis2-template.nc"
SOCA_ANA="${SOCA_TEST}/Data/ocn.3dvar.an.2018-04-15T00:00:00Z.nc"
if ( ${SOCA_DA_IC} ) then
   # Run 3DVAR
   #----------
   #singularity exec ${JCSDA_SING} ${SOCA_BIN}/soca_3dvar.x ${SOCA_INPUT}/3dvar_test.yml

   # Put analysis in restart file
   #-----------------------------
   echo "Making copy of bkg and ana"
   cp ${SOCA_BKG} ${SOCA_TEST}/MOM-ana.res.nc
   cp ${SOCA_ANA} ${SOCA_TEST}/3dvar.an.tmp.nc

   #echo "Rename variable in ana file"
   #ncrename -v temp,Temp 3dvar.an.tmp.nc
   #ncrename -v salt,Salt 3dvar.an.tmp.nc

   echo "Dump analysis variables into bkg restart file"
   ncks -A -v Temp 3dvar.an.tmp.nc MOM-ana.res.nc
   ncks -A -v Salt 3dvar.an.tmp.nc MOM-ana.res.nc

   SCRATCH=${SOCA_TEST}/scratch
else
    # No DA
    #-----------------------------
    echo "Making copy of bkg IC"
    #cp ${SOCA_TEST}/INPUT/MOM.res.nc ${SOCA_TEST}/MOM.res.nc
    cp ${SOCA_BKG} ${SOCA_TEST}/MOM.res.nc    
    SCRATCH=${SOCA_TEST}/scratch_noda
fi

# Setup model
#------------
rm -rf ${SCRATCH}
mkdir  ${SCRATCH}

# Copy resource files in scratch
cp ${SOCA_MODEL_RSC}/* ${SCRATCH} 
mkdir ${SCRATCH}/INPUT
mkdir ${SCRATCH}/RESTART
#cp ${SOCA_MODEL_RSC}/INPUT/ocean_hgrid.nc ${SCRATCH}/INPUT/
#cp ${SOCA_MODEL_RSC}/INPUT/topog.nc ${SCRATCH}/INPUT/
cp ./INPUT/ice_model.res.nc ${SCRATCH}/INPUT/
ln -s ${MOM6_MINSCRATCH}/INPUT/* ${SCRATCH}/INPUT/
cp ${MOM6_MINSCRATCH}/* ${SCRATCH}

# Get restart
rm ${SCRATCH}/INPUT/MOM.res.nc

if ( ${SOCA_DA_IC} ) then
   #if DA on
   cp ${SOCA_TEST}/MOM-ana.res.nc ${SCRATCH}/INPUT/MOM.res.nc
   cp ${SOCA_TEST}/INPUT/ice_model.res.nc ${SCRATCH}/INPUT/ice_model.res.nc   
else
    #if no DA
    cp ${SOCA_TEST}/MOM.res.nc ${SCRATCH}/INPUT/MOM.res.nc
    cp ${SOCA_TEST}/INPUT/ice_model.res.nc ${SCRATCH}/INPUT/ice_model.res.nc
fi

# Start MOM6-SIS2
#----------------
cd ${SCRATCH}
echo $SOCA_BIN   
singularity exec ${JCSDA_SING_old} mpirun -np 4 --oversubscribe /home/gvernier/Sandboxes/MOM6-examples/build/gnu/ice_ocean_SIS2/repro/MOM6
