# (C) COPYRIGHT 2018-2023 UCAR
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

# Use external jedi-cmake or build in bundle
if(DEFINED ENV{jedi_cmake_ROOT})
  include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )
else ()
  ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/jcsda-internal/jedi-cmake.git" UPDATE BRANCH develop RECURSIVE )
  include( jedicmake/cmake/Functions/git_functions.cmake  )
endif()


#===================================================================================================
# other optional repositories
#===================================================================================================
option(BUILD_OASIM "download and build OASIM ocean color forward operator" OFF)
if ( BUILD_OASIM )
  ecbuild_bundle( PROJECT oasim         GIT "https://github.com/JCSDA-internal/oasim.git"           UPDATE BRANCH develop )
endif ()

#===================================================================================================
# required repositories
#===================================================================================================
ecbuild_bundle( PROJECT gsw             GIT "https://github.com/jcsda-internal/GSW-Fortran.git"     UPDATE BRANCH develop )
ecbuild_bundle( PROJECT oops            GIT "https://github.com/jcsda-internal/oops.git"            UPDATE BRANCH develop )
ecbuild_bundle( PROJECT vader           GIT "https://github.com/jcsda-internal/vader.git"           UPDATE BRANCH develop )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/jcsda-internal/saber.git"           UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/jcsda-internal/ioda.git"            UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/jcsda-internal/ufo.git"             UPDATE BRANCH develop )
ecbuild_bundle( PROJECT soca            SOURCE "../" )

#===================================================================================================
# optional test data for the upstream repos
#===================================================================================================
option(TEST_IODA "download and ioda-data repo needed for ioda ctests" OFF)
if ( TEST_IODA )
  find_branch_name(REPO_DIR_NAME ioda)
  if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
    ecbuild_bundle( PROJECT ioda-data GIT "https://github.com/JCSDA-internal/ioda-data.git" BRANCH develop UPDATE )
  endif()
endif ()

option(TEST_UFO "download and ufo-data repo needed for ufo ctests" OFF)
if ( TEST_UFO )
  find_branch_name(REPO_DIR_NAME ufo)
  if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
    ecbuild_bundle( PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" BRANCH develop UPDATE )
  endif()
endif ()

option(TEST_SABER "download and saber-data repo needed for ioda ctests" OFF)
if ( TEST_SABER )
  find_branch_name(REPO_DIR_NAME saber)
  if( NOT DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES} AND NOT DEFINED ${GIT_TAG_FUNC} )
    ecbuild_bundle( PROJECT saber-data GIT "https://github.com/JCSDA-internal/saber-data.git" BRANCH develop UPDATE )
  endif()
endif ()

#===================================================================================================

ecbuild_bundle_finalize()
