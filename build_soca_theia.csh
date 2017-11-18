#!/bin/csh -f                                                                                                                                                               
                                                                                                                                                                            
# Clean modules                                                                                                                                                             
source $MODULESHOME/init/csh                                                                                                                                                
                                                                                                                                                                            
module purge                                                                                                                                                                
                                                                                                                                                                            
# Load Intel compilers, NetCDF and HDF5 libraries                                                                                                                           
module load intel/15.6.233                                                                                                                                                  
module load impi/5.1.2.150                                                                                                                                                  
module load hdf5/1.8.14 netcdf/4.3.0                                                                                                                                        
                                                                                                                                                                            
# Load cmake                                                                                                                                                                
module use -a /contrib/modulefiles                                                                                                                                          
module load cmake                                                                                                                                                           
                                                                                                                                                                            
# Load Boost and Eigen                                                                                                                                                      
module use -a /contrib/da/modulefiles                                                                                                                                       
module load boost                                                                                                                                                           
module load eigen                                                                                                                                                           
                                                                                                                                                                            
module list                                                                                                                                                                 
                                                                                                                                                                            
# Need eigen3 library                                                                                                                                                       
setenv EIGEN3_INCLUDE_DIR $EIGEN_ROOT                                                                                                                                       
                                                                                                                                                                            
# Need NETCDF library                                                                                                                                                       
setenv NETCDF_LIBRARIES "${NETCDF}/lib/libnetcdf.a;${NETCDF}/lib/libnetcdff.a"                                                                                              
                                                                                                                                                                            
# Need FMS library (For FV3GFS model only)                                                                                                                                  
#setenv FMS_ROOT "/scratch4/NCEPDEV/global/save/Rahul.Mahajan/git/GFDL/FMS"                                                                                                 
#setenv FMS_LIBRARIES "${FMS_ROOT}/build/libfms.a"                                                                                                                          
#setenv FMS_INCLUDES  "${FMS_ROOT}/build;${FMS_ROOT}/include"                                                                                                               
                                                                                                                                                                            
                                                                                                                                                                            
# Define source and build directories                                                                                                                                       
                                                                                                                                                                            
setenv SRC_OOPS "/scratch4/NCEPDEV/ocean/scrub/Guillaume.Vernieres/JEDI/TEST/jedi/"                                                                                         
setenv SRC_MODEL "/scratch4/NCEPDEV/ocean/scrub/Guillaume.Vernieres/JEDI/TEST/soca/"                                                                                   
setenv BUILD "/scratch4/NCEPDEV/ocean/scrub/Guillaume.Vernieres/JEDI/TEST/jedi/build/"                                                                                      
setenv PATH ${PATH}:${SRC_OOPS}/ecbuild/bin                                                                                                                                 
#set path = ( ${SRC_OOPS}/ecbuild/bin $path )                                                                                                                               
rm -rf ${BUILD}/soca; mkdir ${BUILD}/soca; cd ${BUILD}/soca                                                                                                  
ecbuild \                                                                                                                                                                   
    -DENABLE_CXX11=ON \                                                                                                                                                     
    -DOOPS_PATH=${BUILD}/oops \                                                                                                                                             
    -DCMAKE_CXX_COMPILER=mpiicpc \                                                                                                                                          
    -DCMAKE_C_COMPILER=mpiicc \                                                                                                                                             
    -DCMAKE_Fortran_COMPILER=mpiifort \                                                                                                                                     
    --build=release \                                                                                                                                                       
    ${SRC_MODEL}                                                                                                                                                            
make -j4                                                                                                                                                                    
                                                                                                                                                                            
#    -DENABLE_CXX11=ON \                                                                                                                                                    
#ecbuild \                                                                                                                                                                  
#    --build=release \                                                                                                                                                      
#    -DCMAKE_CXX_COMPILER=mpiicpc \                                                                                                                                         
#    -DCMAKE_C_COMPILER=mpiicc \                                                                                                                                            
#    -DCMAKE_Fortran_COMPILER=mpiifort \                                                                                                                                    
#    -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_SYSTEM_PATHS=ON \                                                                                                                  
#    ${SRC}                                                                                                                                                                 
                                                                                                                                                                            
# Compile                                                                                                                                                                   
#make VERBOSE=YES -j4                                                                                                                                                       
                                                                                                                                                                            
exit 0                                                                                                                                                                      
                                                                                                                                                                            
                                                                                                                                                                            
                                                                                                                                                                            
     
