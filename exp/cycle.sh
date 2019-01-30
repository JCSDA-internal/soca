#!/bin/sh

# Set date for IC
YYYY=2018
MM=07
DD=22
HH=12

# DA window
WINDOW_LENGTH=PT24H

# Initialize B
INIT_B=false 

# Change dates in DA config file
WINDOW_BEGIN=${YYYY}-${MM}-${DD}T${HH}:00:00Z
cp ./expinput/3dvarfgat.yml 3dvarfgat-${WINDOW_BEGIN}.yml 
sed -i "s/WINDOW_BEGIN/${WINDOW_BEGIN}/g" 3dvarfgat-${WINDOW_BEGIN}.yml 
sed -i "s/WINDOW_LENGTH/${WINDOW_LENGTH}/g" 3dvarfgat-${WINDOW_BEGIN}.yml 

# Update MOM6's IC date
cp input-ymd.nml input.nml
sed -i "s/YYYY/${YYYY}/g" input.nml
sed -i "s/MM/${MM}/g" input.nml 
sed -i "s/DD/${DD}/g" input.nml
sed -i "s/HH/${HH}/g" input.nml

# Prepare obs file
WINDOW_MID='2018-07-23T00:00:00Z'

#Jason-3
OBS_ADT_J3=ADT-J3-${WINDOW_MID}
sed -i "s/OBS_ADT_J3/${OBS_ADT_J3}/g" 3dvarfgat-${WINDOW_BEGIN}.yml
cp /home/gvernier/Data/FNMOC/ioda-converters/src/marine/altimeter/adt-ioda-j3-$YYYY$MM$DD$HH.nc ./Data/${OBS_ADT_J3}.nc

#Cryosat-2
OBS_ADT_C2=ADT-C2-${WINDOW_MID}
sed -i "s/OBS_ADT_C2/${OBS_ADT_C2}/g" 3dvarfgat-${WINDOW_BEGIN}.yml
cp /home/gvernier/Data/FNMOC/ioda-converters/src/marine/altimeter/adt-ioda-c2-$YYYY$MM$DD$HH.nc ./Data/${OBS_ADT_C2}.nc

#Cryosat-2
OBS_ADT_SA=ADT-SA-${WINDOW_MID}
sed -i "s/OBS_ADT_SA/${OBS_ADT_SA}/g" 3dvarfgat-${WINDOW_BEGIN}.yml
cp /home/gvernier/Data/FNMOC/ioda-converters/src/marine/altimeter/adt-ioda-sa-$YYYY$MM$DD$HH.nc ./Data/${OBS_ADT_SA}.nc

# Initialize static B
if ( ${INIT_B} ) then
   mpirun -np 4 --oversubscribe ../../bin/soca_staticbinit.x ./expinput/static_SocaError_init.yml
fi

# Run DA
mpirun -np 4 --oversubscribe ../../bin/soca_3dvar.x 3dvarfgat-${WINDOW_BEGIN}.yml

# Update IC
