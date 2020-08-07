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
  ecbuild_bundle( PROJECT eckit         GIT "https://github.com/JCSDA/eckit.git"           UPDATE TAG 784fe3d )
endif ()

option(BUILD_CRTM  "download and build CRTM")
if ( BUILD_CRTM )
  ecbuild_bundle( PROJECT crtm          GIT "https://github.com/JCSDA/crtm.git"            UPDATE TAG 747fc26 )
endif ()

option(BUILD_METIS "download and build METIS (required for points remapping in BUMP)")
if ( BUILD_METIS )
  ecbuild_bundle( PROJECT metis         GIT "https://github.com/JCSDA/metis.git"           UPDATE BRANCH develop )
endif ()


# required repositories

ecbuild_bundle( PROJECT fckit           GIT "https://github.com/JCSDA/fckit.git"           UPDATE TAG 8502aae )
ecbuild_bundle( PROJECT atlas           GIT "https://github.com/JCSDA/atlas.git"           UPDATE TAG e0b7709 )
ecbuild_bundle( PROJECT oops            GIT "https://github.com/JCSDA/oops.git"            UPDATE TAG f99dda8 )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/JCSDA/saber.git"           UPDATE TAG a78365a )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/JCSDA/ioda.git"            UPDATE TAG b88ab53 )
ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/JCSDA/ioda-converters.git" UPDATE TAG f6e4a1b )
ecbuild_bundle( PROJECT gsw             GIT "https://github.com/JCSDA/GSW-Fortran.git"     UPDATE TAG 8987bf8 )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/JCSDA/ufo.git"             UPDATE TAG 9e9201c )
ecbuild_bundle( PROJECT fms             GIT "https://github.com/JCSDA/FMS.git"             UPDATE TAG 5bbd9e6 )
ecbuild_bundle( PROJECT mom6            GIT "https://github.com/JCSDA/MOM6.git"            UPDATE TAG 4469d6d RECURSIVE)
ecbuild_bundle( PROJECT soca-config     GIT "https://github.com/JCSDA/soca-config.git"     UPDATE TAG 7415936 )
ecbuild_bundle( PROJECT soca            SOURCE "../" )

ecbuild_bundle_finalize()
