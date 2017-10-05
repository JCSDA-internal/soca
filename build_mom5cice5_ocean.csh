#!/bin/bash -f

#export CC=gcc

export SRC=$PWD
export BUILD=${SRC}/build

# Remove anaconda from path, reduce to bare minimum.
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Path to JCSDA common_libs
export COMMON_LIBS=/home/gvernier/Sandboxes/common_libs/

# Need eigen
export EIGEN3_INCLUDE_DIR=/usr/include/eigen3/

# Need boost library
export BOOST_ROOT=/home/gvernier/Sandboxes/boost_1_64_0/

# Need openmpi
export OPENMPI_BINDIR=${COMMON_LIBS}/bin


export SRC_OOPS="/home/gvernier/Sandboxes/JEDI/jedi/"
export SRC_MODEL="/home/gvernier/Sandboxes/JEDI/mom5cice5/"
export BUILD=${SRC_OOPS}"build"
export PATH=${PATH}:${SRC_OOPS}/ecbuild/bin
#set path = ( ${SRC_OOPS}/ecbuild/bin $path )
echo $SRC_OOPS
rm -rf ${BUILD}/mom5cice5; mkdir ${BUILD}/mom5cice5; cd ${BUILD}/mom5cice5
    #-DCMAKE_C_COMPILER=mpiCC \
ecbuild \
    -DENABLE_CXX11=ON \
    -DOOPS_PATH=${BUILD}/oops \
    -DCMAKE_CXX_COMPILER=mpiCC \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    --build=release \
    ${SRC_MODEL}
make VERBOSE=1 -j4                                                                                                                                                                    

#######################
#rm -rf ${BUILD}; mkdir ${BUILD}; cd ${BUILD}
#export NETCDF_LIBRARIES="/lib/libnetcdff.a;/lib/libnetcdf.a"
#export PATH=${PATH}:${SRC}/ecbuild/bin:$EIGEN3_INCLUDE_DIR

# configure
#ecbuild \
#    --build=release \
#    -D_GLIBCXX_USE_CXX11_ABI=0 \
#    -DENABLE_CXX11=ON \
#    -DNETCDF_Fortran=ON \
#    -DCMAKE_CXX_COMPILER=${OPENMPI_BINDIR}/mpicxx \
#    -DCMAKE_C_COMPILER=${OPENMPI_BINDIR}/mpicc \
#    -DCMAKE_Fortran_COMPILER=${OPENMPI_BINDIR}/mpifort \
#    -DBOOST_ROOT=$BOOST_ROOT \
#    -DBoost_NO_SYSTEM_PATHS=ON ${SRC}

# Compile
#make VERBOSE=1 -j4

exit 0
