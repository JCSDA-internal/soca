#!/bin/bash -f
export SRC="/home/gvernier/Sandboxes/codesprint-ufo/jedi/code/"
export SRC_OOPS=${SRC}"jedi-bundle/oops/"
export SRC_MODEL=${SRC}"soca/"
export BUILD="/home/gvernier/Sandboxes/codesprint-ufo/jedi/build"
rm -rf ${BUILD}/soca; mkdir ${BUILD}/soca; cd ${BUILD}/soca
ecbuild -DOOPS_PATH=${BUILD}/oops -DIODA_PATH=${BUILD}/ioda -DUFO_PATH=${BUILD}/ufo ${SRC_MODEL}
make VERBOSE=1 -j1
exit 0
