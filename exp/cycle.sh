#!/bin/sh

# Set date for IC
CDATE=$1  # 2018
WINDOW_LENGTH=PT$2H # DA window


YYYY=$(echo $CDATE | cut -c1-4)
MM=$(echo $CDATE | cut -c5-6)
DD=$(echo $CDATE | cut -c7-8)
HH=$(echo $CDATE | cut -c9-10)

export OMP_NUM_THREADS=1

# Initialize B
INIT_B=false
if ( ${INIT_B} ) then
  mpirun ../../bin/soca_staticbinit.x ./expinput/static_SocaError_init.yml
fi

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

#Jason-3
OBS_ADT_J3=ADT-J3-${CDATE}
sed -i "s/OBS_ADT_J3/${OBS_ADT_J3}/g" 3dvarfgat-${WINDOW_BEGIN}.yml
cp adt-ioda-j3-${CDATE}.nc ./Data/${OBS_ADT_J3}.nc

#Cryosat-2
OBS_ADT_C2=ADT-C2-${CDATE}
sed -i "s/OBS_ADT_C2/${OBS_ADT_C2}/g" 3dvarfgat-${WINDOW_BEGIN}.yml
cp adt-ioda-c2-${CDATE}.nc ./Data/${OBS_ADT_C2}.nc

#Sentinel
OBS_ADT_SA=ADT-SA-${YYYY}${MM}${DD}${MM}${HH}
sed -i "s/OBS_ADT_SA/${OBS_ADT_SA}/g" 3dvarfgat-${WINDOW_BEGIN}.yml
cp adt-ioda-sa-${CDATE}.nc ./Data/${OBS_ADT_SA}.nc

# Run DA
mpirun ../../bin/soca_3dvar.x 3dvarfgat-${WINDOW_BEGIN}.yml

# Update IC
mv ./RESTART/MOM.res.nc ./INPUT/
