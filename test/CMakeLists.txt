list( APPEND soca_test_input
  testinput/3dhyb.yml
  testinput/3dhybfgat.yml
  testinput/3dvar_godas.yml
  testinput/3dvar_soca.yml
  testinput/3dvarfgat.yml
  testinput/checkpointmodel.yml
  testinput/dirac_bump_cov.yml
  testinput/dirac_horizfilt.yml
  testinput/dirac_soca_cov.yml
  testinput/dirac_socahyb_cov.yml
  testinput/enshofx.yml
  testinput/enspert.yml
  testinput/ensvariance.yml
  testinput/forecast_identity.yml
  testinput/forecast_identity_cice.yml
  testinput/forecast_mom6.yml
  testinput/geometry.yml
  testinput/gridgen.yml
  testinput/hofx.yml
  testinput/hofx3d.yml
  testinput/increment.yml
  testinput/linearmodel.yml
  testinput/marine2ioda.yml
  testinput/model.yml
  testinput/parameters_bump_cov_lct.yml
  testinput/parameters_bump_cov_nicas.yml
  testinput/parameters_bump_loc.yml
  testinput/state.yml
  testinput/static_socaerror_init.yml
  testinput/static_socaerror_test.yml
  testinput/varchange_balance.yml
  testinput/varchange_bkgerrfilt.yml
  testinput/varchange_bkgerrgodas.yml
  testinput/varchange_bkgerrsoca.yml
  testinput/varchange_horizfilt.yml
  testinput/varchange_vertconv.yml
)

list( APPEND soca_test_ref
  testref/3dhyb.test
  testref/3dhybfgat.test
  testref/3dvar_godas.test
  testref/3dvar_soca.test
  testref/3dvarfgat.test
  testref/checkpointmodel.test
  testref/dirac_bump_cov.test
  testref/dirac_horizfilt.test
  testref/dirac_soca_cov.test
  testref/dirac_socahyb_cov.test
  testref/enshofx.test
  testref/enspert.test
  testref/ensvariance.test
  testref/forecast_identity_cice.test
  testref/forecast_identity.test
  testref/forecast_mom6.test
  testref/hofx.test
  testref/hofx3d.test
  testref/parameters_bump_cov_lct.test
)

list( APPEND soca_model_param
  Data/72x35x25/diag_table
  Data/72x35x25/field_table
  Data/72x35x25/ice.bkgerror.nc
  Data/72x35x25/input.nml
  Data/72x35x25/MOM_input
  Data/72x35x25/ocn.bkgerror.nc
  Data/godas_sst_bgerr.nc
  Data/rossrad.dat
)

list( APPEND soca_model_restarts
  Data/72x35x25/INPUT/forcing_daily.nc
  Data/72x35x25/INPUT/forcing_monthly.nc
  Data/72x35x25/INPUT/hycom1_25.nc
  Data/72x35x25/INPUT/ice_model.res.nc
  Data/72x35x25/INPUT/cice.res.nc
  Data/72x35x25/INPUT/layer_coord25.nc
  Data/72x35x25/INPUT/MOM.res.nc
  Data/72x35x25/INPUT/sfc.res.nc
  Data/72x35x25/INPUT/ocean_hgrid.nc
  Data/72x35x25/INPUT/ocean_topog.nc
)

list( APPEND soca_obs
  Data/Obs/adt.nc
  Data/Obs/icec.nc
  Data/Obs/prof.nc
  Data/Obs/sss.nc
  Data/Obs/sst.nc
)

# Create Data directory for test input/output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
foreach(FILENAME ${soca_test_input} ${soca_test_ref})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for model input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/INPUT)
foreach(FILENAME ${soca_model_restarts})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/INPUT/${filename} )
endforeach(FILENAME)

# MOM's resource files
foreach(FILENAME ${soca_model_param})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)


# MOM's output directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RESTART)

ecbuild_add_resources( TARGET   soca_test_scripts
                       SOURCES_PACK
                       ${soca_test_input}
                       )

# Link to SOCA obs
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
foreach(FILENAME ${soca_obs})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )
endforeach(FILENAME)

# Create directory for BUMP output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)


#####################################################################

# Configurable parameters for the subsequent tests

set( SOCA_TESTS_MAX_MPI 2 CACHE STRING
  "Maximum number of PEs to use for SOCA's MPI based tests.\
  NOTE: should leave this at 2 to ensure all tests pass.")

set( SOCA_TESTS_DEFAULT_TOL "1.0e-12;0" CACHE STRING
  "Default test tolerances for compare.py based tests. (float_tol;int_tol) ")

set( SOCA_TESTS_FORC_DEFAULT_TOL OFF CACHE BOOL
  "If true, the specified test tolerances are ignored and\
  SOCA_TESTS_DEFAULT_TOL is used" )

set( SOCA_TESTS_FORC_TRAPFPE OFF CACHE BOOL
  "If true, the specified per-test trapping on/off flag is ignored\
  and is forced to ON for all tests" )


#-------------------------------------------------------------------------------
# The following is a wrapper to simplify the generation of tests.
# There are two types of tests:
#  1) soca executables (EXE must be given)
#  2) interface tests  (SRC must be given)
# Soca executables by default will also run a second step where the output test
# log stream results are compared with reference values. The interface tests
# do not do this.
#
# Arguments:
#  NOCOMPARE - if flag is present, no compare test is performed
#  NAME      - the name of the test (after test_soca is prepended), as well
#              as the base name of the input confg file
#  EXE       - name of executable to be used.
#  SRC       - The name of the unit test source file to compile
#  CFG       - The name of the yaml file, if the default testinput/${NAME}.yaml
#              is not to be used
#  TOL       - The tolerances to use for the test "float_tolerance int_tolerance"
#              If none are given, ${SOCA_TESTS_DEFAULT_TOL} is used. Only used
#              if an EXE with no NOCOMPARE flag set.
# TEST_DEPENDS - list of tests this test depends on
#-------------------------------------------------------------------------------
function(soca_add_test)
  # parse the passed arguments
  set(prefix     ARG)
  set(novals     NOCOMPARE NOTRAPFPE)
  set(singlevals NAME EXE SRC CFG)
  set(multivals  TEST_DEPENDS TOL)
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

  # set the tolerances to defaults if none given, or if we want to override
  IF(NOT ARG_TOL OR SOCA_TESTS_FORC_DEFAULT_TOL)
    SET(ARG_TOL ${SOCA_TESTS_DEFAULT_TOL})
  ENDIF()

  # determine if floating point error trapping should be set
  if ( ARG_NOTRAPFPE AND NOT SOCA_TESTS_FORC_TRAPFPE)
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=0")
  else()
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=1")
  endif()

  # determine the default config file name
  if ( ARG_CFG )
    set ( CONFIG_FILE testinput/${ARG_CFG} )
  else()
    set ( CONFIG_FILE testinput/${ARG_NAME}.yml )
  endif()

  # Are we building a unit test / or running a soca executable?
  if ( ARG_SRC )
    # building a unit test, therfore also assume no compare step
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      SOURCES executables/${ARG_SRC}
                      ARGS    ${CONFIG_FILE}
                      LIBS    soca
                      MPI     ${SOCA_TESTS_MAX_MPI}
                      ENVIRONMENT ${TRAPFPE_ENV}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
    set ( ARG_NOCOMPARE TRUE)

  else()
    # running a soca executable
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      TYPE    EXE
                      COMMAND "${CMAKE_BINARY_DIR}/bin/${ARG_EXE}"
                      ARGS    ${CONFIG_FILE}
                              testoutput/${ARG_NAME}.log
                      MPI     ${SOCA_TESTS_MAX_MPI}
                      ENVIRONMENT ${TRAPFPE_ENV}
                      DEPENDS ${ARG_EXE}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
  endif()

  # compare the answers, if desired
  IF(NOT ARG_NOCOMPARE)
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}_compare
                      TYPE    SCRIPT
                      COMMAND "${CMAKE_BINARY_DIR}/bin/compare.py"
                      ARGS    testoutput/${ARG_NAME}.log
                              testref/${ARG_NAME}.test
                              ${ARG_TOL}
                     TEST_DEPENDS test_soca_${ARG_NAME})
  ENDIF()
endfunction()



#-------------------------------------------------------------------------------
# Tests that create data other tests will use
#-------------------------------------------------------------------------------

# Test of grid generation
soca_add_test( NAME gridgen
               EXE  soca_gridgen.x
               NOCOMPARE )

# Test of read/write restart for MOM6&SIS2, model advance = Id
soca_add_test( NAME forecast_identity
               EXE  soca_forecast.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen)

# Test of read/write restart for MOM6&CICE, model advance = Id
soca_add_test( NAME forecast_identity_cice
               EXE  soca_forecast.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Test mom6 forecast
soca_add_test( NAME forecast_mom6
               EXE  soca_forecast.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Initialize static error covariance model defined in soca, used for 3DVAR
soca_add_test( NAME static_socaerror_init
               EXE  soca_staticbinit.x
               NOCOMPARE
               TEST_DEPENDS test_soca_gridgen )

# ensemble generation
soca_add_test( NAME enspert
               EXE  soca_enspert.x
               TOL 7.0e-4 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_static_socaerror_init)

# TODO Test below can take a few hours to complete.
# Revive once we have tiered testing?
#soca_add_test( NAME covar_lengthscale
#               EXE  soca_staticbinit.x
#               TEST_DEPENDS test_soca_gridgen)


# Tests of class interfaces
#-------------------------------------------------------------------------------
soca_add_test( NAME geometry
               SRC  TestGeometry.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME state
               SRC  TestState.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME model
               SRC  TestModel.cc
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME modelaux
               SRC  TestModelAuxControl.cc
               CFG  model.yml
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME increment
               SRC  TestIncrement.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME errorcovariance
               SRC  TestErrorCovariance.cc
               CFG  static_socaerror_test.yml
               TEST_DEPENDS test_soca_gridgen
                            test_soca_static_socaerror_init )

soca_add_test( NAME linearmodel
               SRC  TestLinearModel.cc
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_static_socaerror_init )

soca_add_test( NAME varchange_balance
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrfilt
               SRC TestVariableChange.cc
                         TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_horizfilt
               SRC TestVariableChange.cc
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrsoca
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrgodas
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_vertconv
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )


#-------------------------------------------------------------------------------
# Tests using the main SOCA applications
#  TEST output stream is compared against existing testoutput/ files
#-------------------------------------------------------------------------------
soca_add_test( NAME ensvariance
               EXE  soca_ensvariance.x
               TOL  3.0e-4 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_enspert )

soca_add_test( NAME parameters_bump_loc
               EXE  soca_parameters.x
               NOCOMPARE
               TEST_DEPENDS test_soca_enspert )

# TODO the next couple of tests are a little slow.
# Revive them once we have tiered testing?
# soca_add_test( NAME parameters_bump_cov_lct
#                EXE  soca_parameters.x
#                TEST_DEPENDS test_soca_enspert )
#
# soca_add_test( NAME parameters_bump_cov_nicas
#                EXE  soca_parameters.x
#                TEST_DEPENDS test_soca_parameters_bump_cov_lct )

# soca_add_test( NAME dirac_bump_cov
#                EXE  soca_dirac.x )


# Dirac tests

soca_add_test( NAME dirac_soca_cov
               EXE  soca_dirac.x
               TOL  8.0e-3 0
               TEST_DEPENDS test_soca_static_socaerror_init)

soca_add_test( NAME dirac_socahyb_cov
               EXE  soca_dirac.x
               TOL  5.0e-3 0
               TEST_DEPENDS test_soca_parameters_bump_loc)

soca_add_test( NAME dirac_horizfilt
               EXE  soca_dirac.x
               NOTRAPFPE
               TOL  5.0e-3 0
               TEST_DEPENDS test_soca_static_socaerror_init)

#ecbuild_add_test( TARGET  test_soca_localization
#                  SOURCES executables/TestLocalization.cc
#                  ARGS    "testinput/localization_test.yml"
#                  MPI     1
#                  LIBS    soca )

#ecbuild_add_test( TARGET test_soca_obserrorcov
#                  SOURCES executables/TestObsErrorCovariance.cc
#                  ARGS "testinput/obserrorcovariance_test.yml"
#                  LIBS    soca )
#

# h(x) executables
soca_add_test( NAME hofx3d
               EXE  soca_hofx3d.x
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME hofx
               EXE  soca_hofx.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# TODO enshofx currently breaks compare.py
# Re-enable the comparison once that's fixed
soca_add_test( NAME enshofx
               EXE  soca_enshofx.x
               NOCOMPARE
               NOTRAPFPE
               TEST_DEPENDS test_soca_enspert )


# variational methods
soca_add_test( NAME 3dvar_soca
               EXE  soca_3dvar.x
               TOL  2.3e-4 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dvar_godas
               EXE  soca_3dvar.x
               TOL  2e-2 0
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dvarfgat
               EXE  soca_3dvar.x
               TOL  2e-4 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dhyb
               EXE  soca_3dvar.x
               TOL  5e-4 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init
                            test_soca_parameters_bump_loc )

soca_add_test( NAME 3dhybfgat
               EXE  soca_3dvar.x
               TOL  4e-3 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init
                            test_soca_parameters_bump_loc )

# soca_add_test( NAME 4dhyb
#                EXE  soca_3dvar.x
#                TEST_DEPENDS test_soca_static_socaerror_init
#                             test_soca_parameters_bump_loc )


# restart checkpointing
soca_add_test( NAME checkpointmodel
               EXE  soca_checkpoint_model.x
               TOL  2e-6 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_3dvar_godas)
