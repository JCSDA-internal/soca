# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


################################################################################
# SOCA
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( soca VERSION 1.0.0 LANGUAGES C CXX Fortran)

find_package(ecbuild 3.3.2 REQUIRED)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( soca_compiler_flags )
################################################################################
# Dependencies
################################################################################
find_package( NetCDF REQUIRED COMPONENTS Fortran)
find_package( gsl-lite REQUIRED HINTS $ENV{gsl_list_DIR})
find_package( eckit  1.11.6   REQUIRED)
find_package( fckit  0.7.0    REQUIRED)
find_package( atlas  0.20.2   REQUIRED)
find_package( oops   1.0.0    REQUIRED)
find_package( saber  1.0.0    REQUIRED)
find_package( ioda   1.0.0    REQUIRED)
find_package( ufo    1.0.0    REQUIRED)
#find_package( fms    2020.4.0 REQUIRED)
#find_package( mom6   2020.4.0 REQUIRED)
find_package(ESMF MODULE REQUIRED STATIC)
find_package( crtm   2.2.3)
  find_package(bacio REQUIRED )
  find_package(sigio REQUIRED )
  find_package(sp REQUIRED )
  find_package(w3emc REQUIRED )
  find_package(w3nco REQUIRED )
  find_package(nemsio REQUIRED )
include_directories( ${NETCDF_INCLUDE_DIRS} ${DEPEND_LIB_ROOT}/include_r8 ${DEPEND_LIB_ROOT}/MOM6/mod  )
include_directories( ${DEPEND_LIB_ROOT}/../ufs-weather-model/src/ufs-weather-model-build/MOM6-interface/mod_solo)


################################################################################
# Sources
################################################################################
set( SOCA_LINKER_LANGUAGE CXX)

if(${UFS_APP} MATCHES "NG-GODAS")
  set( UFS_LIBS 
      soca
      cdeps
      cmeps
      piof 
      pioc
      stochastic_physics
      cice
      esmf
  )
elseif(${UFS_APP} MATCHES "S2S")
  set( UFS_LIBS 
      soca
      cmeps
      piof 
      pioc
      fv3atm
      fv3
      fv3ccpp
      ccpp_framework
      ccpp_physics
      stochastic_physics
      FMS::fms_r8
      cice
      esmf
  )
else()
    message(FATAL_ERROR "soca UFS_LIBS not configured for UFS_APP ${UFS_APP}")
endif()

add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( tutorial )

# DH* 20221120 - temporary until debugging the build is done - use relative dir, also mkdir ...
include(FetchContent)

FetchContent_Declare(
  testdata
# URL      https://ftp.emc.ncep.noaa.gov/jcsda/WDQMS/NCEP/UFS-JEDI-DATA/test.tar.gz
# URL_HASH MD5=822ffcd68e5d35541fa183d6811d690b
  URL      https://ftp.emc.ncep.noaa.gov/jcsda/WDQMS/NCEP/UFS-JEDI-DATA/ufs.tar.gz
  URL_HASH MD5=762e54986481c01e3a939edbe57b0ec4
)
FetchContent_Populate(
  testdata
  #SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../cache/soca_test_data_762e54986481c01e3a939edbe57b0ec4
  SOURCE_DIR /Users/heinzell/work/ufs-bundle/20221114/ufs-bundle/cache/soca_test_data_762e54986481c01e3a939edbe57b0ec4
)
FetchContent_MakeAvailable(testdata)
# *DH 20221120

# Build Doxygen documentation
add_subdirectory( docs )

if(ECBUILD_INSTALL_FORTRAN_MODULES)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR} )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME soca )

# print the summary of the configuration
ecbuild_print_summary()
