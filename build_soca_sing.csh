#!/bin/bash -f
export SRC="/home/gvernier/Sandboxes/codesprint-ufo/jedi/code/"
export SRC_OOPS=${SRC}"jedi-bundle/oops/"
export SRC_MODEL=${SRC}"soca/"
export BUILD="/home/gvernier/Sandboxes/codesprint-ufo/jedi/build"


export MOM6_SRC="/home/gvernier/Sandboxes/MOM6-examples/src/"
export MOM6_PATH="/home/gvernier/Sandboxes/MOM6-examples/build/gnu/ice_ocean_SIS2/repro/"
export MOM6_LIBRARIES="${MOM6_PATH}libmom6.a"
export MOM6_INCLUDE_DIRS=$MOM6_PATH
export MOM6_HEADER_DIRS=${MOM6_SRC}/MOM6/config_src/dynamic/
export MOM6_MACROS_DIRS=${MOM6_SRC}/MOM6/src/framework/

export FMS_PATH="/home/gvernier/Sandboxes/MOM6-examples/build/gnu/shared/repro/"
export FMS_LIBRARIES="$FMS_PATH/libfms.a"
export FMS_INCLUDE_DIRS="$FMS_PATH"

export MPI_INCLUDE_DIRS="/usr/local/lib/"
export MPI_LIBRARIES="/usr/local/lib/libmpi.so"

rm -rf ${BUILD}/soca; mkdir ${BUILD}/soca; cd ${BUILD}/soca
ecbuild -DOOPS_PATH=${BUILD}/oops -DIODA_PATH=${BUILD}/ioda -DUFO_PATH=${BUILD}/ufo -DMOM6_LIBRARIES=${MOM6_LIBRARIES} -DMOM6_INCLUDE_DIRS=${MOM6_INCLUDE_DIRS} -DMOM6_HEADER_DIRS=${MOM6_HEADER_DIRS} -DMOM6_MACROS_DIRS=${MOM6_MACROS_DIRS} -DFMS_LIBRARIES=${FMS_LIBRARIES} -DFMS_INCLUDE_DIRS=${FMS_INCLUDE_DIRS} -DMPI_LIBRARIES=${MPI_LIBRARIES} -DMPI_INCLUDE_DIRS=${MPI_INCLUDE_DIRS} ${SRC_MODEL}
make VERBOSE=1 -j1

exit 0
