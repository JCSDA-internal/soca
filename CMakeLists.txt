# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})
#set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake)
include( ecbuild_bundle )
project( soca C CXX Fortran)

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ECBUILD_DEFAULT_BUILD_TYPE Release )
include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )
################################################################################
# Project
################################################################################

ecbuild_declare_project()
ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()
set( SOCA_LINKER_LANGUAGE CXX )

################################################################################
# Dependencies
################################################################################
# Boost
set( Boost_MINIMUM_VERSION "1.47" )
find_package( Boost ${Boost_MINIMUM_VERSION} REQUIRED
COMPONENTS unit_test_framework )
include_directories( ${Boost_INCLUDE_DIR} )

# ioda
ecbuild_use_package( PROJECT ioda VERSION 0.1.0 REQUIRED )
include_directories( ${IODA_INCLUDE_DIRS} )

# ufo
ecbuild_use_package( PROJECT ufo VERSION 0.1.0 REQUIRED )
include_directories( ${UFO_INCLUDE_DIRS} )
MESSAGE( STATUS "ufo include dirs: " ${UFO_INCLUDE_DIRS} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.1.0 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )
MESSAGE( STATUS "oops include dirs: " ${OOPS_INCLUDE_DIRS} )

# fms
ecbuild_use_package( PROJECT fms REQUIRED )
include_directories( ${FMS_INCLUDE_DIRS} )
MESSAGE( STATUS "fms include dirs: " ${FMS_INCLUDE_DIRS} )

# mom6
ecbuild_use_package( PROJECT mom6 REQUIRED )
include_directories( ${MOM6_INCLUDE_DIRS} )
MESSAGE( STATUS "mom6 include dirs: " ${MOM6_INCLUDE_DIRS} )

################################################################################
# Export package info
################################################################################
list( APPEND SOCA_TPLS ufo ioda oops fms mom6)
foreach( _tpl ${SOCA_TPLS} )
string( TOUPPER ${_tpl} TPL )
list( APPEND SOCA_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
list( APPEND SOCA_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
list( APPEND SOCA_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################
include( soca_compiler_flags )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_subdirectory( mains )
add_subdirectory( src )
add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME soca )

# print the summary of the configuration
ecbuild_print_summary()

ecbuild_add_test( TARGET soca_coding_norms
                  TYPE SCRIPT
		  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
		  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

#cmake_minimum_required(VERSION 3.6)

#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )

#add_subdirectory( mains )
#add_subdirectory( src )
#add_subdirectory( test )

################################################################################
# Build documentation
################################################################################

option(BUILD_DOC "Build documentation" ON)
find_package(Doxygen)
if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")
    add_custom_target( doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
	VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
