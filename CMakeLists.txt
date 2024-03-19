# (C) Copyright 2017-2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


################################################################################
# SOCA
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( soca VERSION 1.7.0 LANGUAGES C CXX Fortran)

find_package(ecbuild 3.3.2 REQUIRED)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( soca_compiler_flags )

################################################################################
# Dependencies
################################################################################
find_package( NetCDF REQUIRED COMPONENTS Fortran C)
find_package( gsl-lite REQUIRED HINTS $ENV{gsl_list_DIR})
find_package( eckit  1.23.0   REQUIRED)
find_package( fckit  0.10.1   REQUIRED)
find_package( atlas  0.35.0   REQUIRED)
find_package( oops   1.9.0    REQUIRED)
find_package( vader  1.6.0    REQUIRED)
find_package( saber  1.9.0    REQUIRED)
find_package( ioda   2.8.0    REQUIRED)
find_package( ufo    1.9.0    REQUIRED)
find_package( fms    2020.4.0 REQUIRED) # NOTE: remove this when we update to latest MOM6
find_package( mom6   2020.4.0 REQUIRED) # NOTE: remove this when we update to latest MOM6

# Optionally look for icepack
find_package( icepack 1.2.0 QUIET)

include_directories( ${NETCDF_INCLUDE_DIRS} )

################################################################################
# Model
# TODO update to allow switching between mom6solo (here) and UFS (external)
################################################################################
# # NOTE enable these when we're ready to update to latest MOM6
# find_package( fms    2023.3.0 REQUIRED COMPONENTS R4 R8)
# add_library(  fms ALIAS FMS::fms_r8)
# add_subdirectory(mom6)

################################################################################
# Sources
################################################################################
set( SOCA_LINKER_LANGUAGE CXX)

add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( tutorial )

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
