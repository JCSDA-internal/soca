# (C) Copyright 2017 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake)
include( ecbuild_bundle )
project( mom5cice5 C CXX Fortran)

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ECBUILD_DEFAULT_BUILD_TYPE Release )
#set( ECBUILD_DEFAULT_BUILD_TYPE Debug )
include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )
################################################################################
# Project
################################################################################

ecbuild_declare_project()
ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()
set( MOM5CICE5_LINKER_LANGUAGE CXX )

################################################################################
# Dependencies
################################################################################
# Boost
set( Boost_MINIMUM_VERSION "1.47" )
find_package( Boost ${Boost_MINIMUM_VERSION} REQUIRED
COMPONENTS unit_test_framework )
include_directories( ${Boost_INCLUDE_DIR} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.1.0 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

################################################################################
# Export package info
################################################################################
list( APPEND MOM5CICE5_TPLS oops )
foreach( _tpl ${MOM5CICE5_TPLS} )
string( TOUPPER ${_tpl} TPL )
list( APPEND MOM5CICE5_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
list( APPEND MOM5CICE5_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
list( APPEND MOM5CICE5_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/model )
add_subdirectory( mains )
add_subdirectory( model )
add_subdirectory( test )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME mom5cice5 )

# print the summary of the configuration
ecbuild_print_summary()

#cmake_minimum_required(VERSION 3.6)

#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/model )

#add_subdirectory( mains )
#add_subdirectory( model )
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
