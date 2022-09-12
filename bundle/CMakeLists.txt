# (C) COPYRIGHT 2018-2020 UCAR
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

#
# bundle for SOCA (ocean only)
#
cmake_minimum_required( VERSION 3.12 )
project( soca-bundle VERSION 0.0.0 LANGUAGES C CXX Fortran )

find_package( ecbuild 3.3.2 REQUIRED)
include( ecbuild_bundle )

ecbuild_bundle_initialize()

ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda-internal/jedi-cmake.git" UPDATE TAG b7c8fb7 RECURSIVE )
include( jedicmake/cmake/Functions/git_functions.cmake  )

# ECMWF repos that you probably dont need to build yourself because they
#  should be in jedi-stack and the containers

option(BUILD_ECKIT "download and build eckit (not needed if in a jedi container)")
if ( BUILD_ECKIT )
  ecbuild_bundle( PROJECT eckit  GIT "https://github.com/ecmwf/eckit.git"  UPDATE TAG 1.16.0 )
endif ()

option(BUILD_FCKIT "download and build fckit (not needed if in a jedi container)")
if (BUILD_FCKIT)
  ecbuild_bundle( PROJECT fckit  GIT "https://github.com/jcsda-internal/fckit.git"  UPDATE BRANCH release-stable )
  # TODO: replace with https://github.com/ecmwf/fckit.git TAG 0.9.2
endif()

option(BUILD_ATLAS "download and build atlas")
if (BUILD_ATLAS)
  ecbuild_bundle( PROJECT atlas  GIT "https://github.com/ecmwf/atlas.git" TAG 0.29.0 )
endif()


# other optional repositories

option(BUILD_FMS "download and build fms" ON)
if ( BUILD_FMS )
  ecbuild_bundle( PROJECT fms           GIT "https://github.com/jcsda/FMS.git"                     UPDATE TAG 1f73914 )
endif ()

option(BUILD_UFSM6C6 "download and install the UFS" OFF)

# required repositories

ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda-internal/oops.git"            UPDATE TAG 5138524 )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda-internal/saber.git"           UPDATE TAG c135fa9 )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda-internal/ioda.git"            UPDATE TAG 99d796f )
ecbuild_bundle( PROJECT gsw             GIT "https://github.com/jcsda-internal/GSW-Fortran.git"     UPDATE TAG a9a0989 )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda-internal/ufo.git"             UPDATE TAG fa3acfd )
ecbuild_bundle( PROJECT mom6            GIT "https://github.com/jcsda-internal/MOM6.git"            UPDATE TAG 51ec489 RECURSIVE )
ecbuild_bundle( PROJECT soca            SOURCE "../" )

ecbuild_bundle_finalize()
