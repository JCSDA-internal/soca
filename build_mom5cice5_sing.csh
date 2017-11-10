#!/bin/bash -f
export SRC="/home/gvernier/Sandboxes/codesprint-ufo/jedi/code/"
export SRC_OOPS=${SRC}"oops-nicas/"
export SRC_MODEL=${SRC}"mom5cice5/"
export BUILD="/home/gvernier/Sandboxes/codesprint-ufo/jedi/build"
#export PATH=${PATH}:${SRC_OOPS}/ecbuild/bin
#set path = ( ${SRC_OOPS}/ecbuild/bin $path )
echo $SRC_OOPS
rm -rf ${BUILD}/mom5cice5; mkdir ${BUILD}/mom5cice5; cd ${BUILD}/mom5cice5
    #-DCMAKE_C_COMPILER=mpiCC \
ecbuild \
    -DENABLE_CXX11=ON \
    -DOOPS_PATH=${BUILD}/oops-nicas \
    -DCMAKE_CXX_COMPILER=mpiCC \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    --build=release \
    ${SRC_MODEL}
make VERBOSE=1 -j4
exit 0
