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

ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda-internal/jedi-cmake.git" UPDATE BRANCH develop )
include( jedicmake/cmake/Functions/git_functions.cmake  )

# ECMWF repos that you probably dont need to build yourself because they
#  should be in jedi-stack and the containers

ecbuild_bundle( PROJECT eckit GIT "https://github.com/ecmwf/eckit.git" TAG 1.16.0 )
ecbuild_bundle( PROJECT fckit GIT "https://github.com/ecmwf/fckit.git" TAG 0.9.2 )
ecbuild_bundle( PROJECT atlas GIT "https://github.com/ecmwf/atlas.git" TAG 0.24.1 )

# other optional repositories

option(BUILD_CRTM  "download and build CRTM")
if ( BUILD_CRTM )
  ecbuild_bundle( PROJECT crtm          GIT "https://github.com/jcsda-internal/crtm.git"            UPDATE BRANCH release/crtm_jedi )
endif ()

option(BUILD_FMS "download and build fms" ON)
if ( BUILD_FMS )
  ecbuild_bundle( PROJECT fms           GIT "https://github.com/jcsda/FMS.git"                     UPDATE BRANCH release-stable )
endif ()


# required repositories

ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda-internal/oops.git"            UPDATE BRANCH feature/loc_refactor )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda-internal/saber.git"           UPDATE BRANCH feature/loc_refactor )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda-internal/ioda.git"            UPDATE BRANCH develop )
ecbuild_bundle( PROJECT gsw             GIT "https://github.com/jcsda-internal/GSW-Fortran.git"     UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda-internal/ufo.git"             UPDATE BRANCH develop )
ecbuild_bundle( PROJECT mom6            GIT "https://github.com/jcsda-internal/MOM6.git"            UPDATE BRANCH main-ecbuild RECURSIVE )
ecbuild_bundle( PROJECT soca            SOURCE "../" )

ecbuild_bundle_finalize()
