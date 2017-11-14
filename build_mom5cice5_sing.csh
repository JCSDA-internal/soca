#!/bin/bash -f
export SRC="/home/gvernier/Sandboxes/codesprint-ufo/jedi/code/"
export SRC_OOPS=${SRC}"oops-nicas/"
export SRC_MODEL=${SRC}"mom5cice5/"
export BUILD="/home/gvernier/Sandboxes/codesprint-ufo/jedi/build"

#IODA
setenv IODA_PATH ${BUILD}/ioda
#export PATH=${PATH}:${SRC_OOPS}/ecbuild/bin
#set path = ( ${SRC_OOPS}/ecbuild/bin $path )
echo $SRC_OOPS
rm -rf ${BUILD}/mom5cice5; mkdir ${BUILD}/mom5cice5; cd ${BUILD}/mom5cice5
#-DCMAKE_C_COMPILER=mpiCC \
#    -DCMAKE_Fortran_COMPILER=mpif90 \
#   -DCMAKE_CXX_COMPILER=mpiCC \
#    -DENABLE_CXX11=ON \
ecbuild \
     -DOOPS_PATH=${BUILD}/oops-nicas \
     -DIODA_PATH=${BUILD}/ioda \
    --build=release \
    ${SRC_MODEL}
#make VERBOSE=1 -j4
exit 0
