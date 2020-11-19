# (C) COPYRIGHT 2018-2020 UCAR
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

#
# bundle for SOCA (ocean only)
#

project( soca-bundle C CXX Fortran )

cmake_minimum_required( VERSION 3.8 FATAL_ERROR )

include( ecbuild_bundle )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()
ecbuild_requires_macro_version( 2.7 )


# Optional repositories

option(BUILD_ECKIT "download and build eckit (not needed if in a jedi container)")
if ( BUILD_ECKIT )
  ecbuild_bundle( PROJECT eckit         GIT "https://github.com/jcsda-internal/eckit.git"           UPDATE TAG b3adb28 )
endif ()

option(BUILD_CRTM  "download and build CRTM")
if ( BUILD_CRTM )
  ecbuild_bundle( PROJECT crtm          GIT "https://github.com/jcsda-internal/crtm.git"            UPDATE TAG 139d4f9 )
endif ()

option(BUILD_FMS "download and build fms" ON)
if ( BUILD_FMS )
  ecbuild_bundle( PROJECT fms           GIT "https://github.com/jcsda-internal/FMS.git"             UPDATE TAG f28abfe )
endif ()

# required repositories

ecbuild_bundle( PROJECT fckit           GIT "https://github.com/jcsda-internal/fckit.git"           UPDATE TAG 8502aae )
ecbuild_bundle( PROJECT atlas           GIT "https://github.com/jcsda-internal/atlas.git"           UPDATE TAG e0b7709 )
ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda-internal/oops.git"            UPDATE TAG 0ff544d )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda-internal/saber.git"           UPDATE TAG 9a2025d )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda-internal/ioda.git"            UPDATE TAG 32bde24 )
ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/jcsda-internal/ioda-converters.git" UPDATE TAG 015160b )
ecbuild_bundle( PROJECT gsw             GIT "https://github.com/jcsda-internal/GSW-Fortran.git"     UPDATE TAG 8987bf8 )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda-internal/ufo.git"             UPDATE TAG f7c599d )
ecbuild_bundle( PROJECT mom6            GIT "https://github.com/jcsda-internal/MOM6.git"            UPDATE TAG a04fdb5 RECURSIVE)
ecbuild_bundle( PROJECT soca            SOURCE "../" )

ecbuild_bundle_finalize()
