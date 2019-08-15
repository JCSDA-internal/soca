list( APPEND soca_test_input
  testinput/3dhyb.yml
  testinput/3dhybfgat.yml
  testinput/3dvar_godas.yml
  testinput/3dvar_soca.yml
  testinput/3dvarfgat.yml
  testinput/balance.yml
  testinput/bkgerrfilt.yml
  testinput/bkgerrgodas.yml
  testinput/bkgerrsoca.yml
  testinput/checkpointmodel.yml
  testinput/dirac_bump_cov.yml
  testinput/dirac_soca_cov.yml
  testinput/dirac_socahyb_cov.yml
  testinput/enshofx.yml
  testinput/enspert.yml
  testinput/ensvariance.yml
  testinput/forecast_identity_cice.yml
  testinput/forecast_identity.yml
  testinput/forecast_mom6.yml
  testinput/geometry_test.yml
  testinput/gridgen.yml
  testinput/hofx.yml
  testinput/hofx3d.yml
  testinput/increment_test.yml
  testinput/linearmodel_test.yml
  testinput/marine2ioda.yml
  testinput/model_test.yml
  testinput/parameters_bump_cov_lct.yml
  testinput/parameters_bump_cov_nicas.yml
  testinput/parameters_bump_loc.yml
  testinput/state_test.yml
  testinput/static_error_init.yml
  testinput/static_error_test.yml
  testinput/vertconv.yml
)

list( APPEND soca_test_ref
  testref/3dhyb.test
  testref/3dhybfgat.test
  testref/3dvar_godas.test
  testref/3dvar_soca.test
  testref/3dvarfgat.test
  testref/checkpointmodel.test
  testref/dirac_bump_cov.test
  testref/dirac_soca_cov.test
  testref/dirac_socahyb_cov.test
  testref/enshofx.test
  testref/enspert.test
  testref/ensvariance.test
  testref/forecast_identity_cice.test
  testref/forecast_identity.test
  testref/forecast_mom6.test
  testref/gridgen.test
  testref/hofx.test
  testref/hofx3d.test
  testref/parameters_bump_lct_nicas.test
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

# Link to IODA/SOCA Marine obs
#file(GLOB ioda_obs "${ioda_SOURCE_DIR}/test/testinput/marine/*.nc*")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
#foreach(FILENAME ${ioda_obs} ${soca_obs})
#     get_filename_component(filename ${FILENAME} NAME )
#     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
#           ${CMAKE_${FILENAME}
#           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )
#endforeach(FILENAME)

foreach(FILENAME ${soca_obs})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )
endforeach(FILENAME)

# Create directory for BUMP output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)


#####################################################################

# The number of PEs to use for MPI enabled tests
set(MPI_SOCA_TESTS 2)

# The following is a wrapper to simplify the generation of
# tests that need to 1) run an mpi executable then 2)
# compare the test log stream results with reference values.
function(soca_add_test)
  set(prefix     ARG)
  set(novals     NOCOMPARE)
  set(singlevals NAME EXE)
  set(multivals  TEST_DEPENDS)
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

  ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                    TYPE    EXE
                    COMMAND "${CMAKE_BINARY_DIR}/bin/${ARG_EXE}"
                    ARGS    testinput/${ARG_NAME}.yml
                            testoutput/${ARG_NAME}.log
                    MPI     ${MPI_SOCA_TESTS}
                    DEPENDS ${ARG_EXE}
                    TEST_DEPENDS ${TEST_DEPENDS})

  IF(NOT ${ARG_NOCOMPARE})
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}_compare
                      TYPE    SCRIPT
                      COMMAND "${CMAKE_BINARY_DIR}/bin/compare.py"
                      ARGS    testoutput/${ARG_NAME}.log
                              testref/${ARG_NAME}.test
                              0.0 0
                     TEST_DEPENDS test_soca_${ARG_NAME})
  endif()
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
               TEST_DEPENDS gridgen)

# Test of read/write restart for MOM6&CICE, model advance = Id
soca_add_test( NAME forecast_identity_cice
               EXE  soca_forecast.x
               TEST_DEPENDS gridgen )

# Test mom6 forecast
soca_add_test( NAME forecast_mom6
               EXE  soca_forecast.x
               TEST_DEPENDS gridgen )

# Initialize static error covariance model defined in soca, used for 3DVAR
soca_add_test( NAME static_error_init
               EXE  soca_staticbinit.x
               NOCOMPARE)

# ensemble generation
soca_add_test( NAME enspert
               EXE  soca_enspert.x)

# #Test below can take a few hours to complete. Commented out for now
# ecbuild_add_test( TARGET test_soca_covar_lengthscale
#                   TYPE SCRIPT
#                     COMMAND ./compare.sh
#                   ARGS
#                       "${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x
#                       testinput/static_scale_binit.yml"
#                       testoutput/covar_lengthscale.test
#                   DEPENDS soca_staticbinit.x)



# Tests of class interfaces
#-------------------------------------------------------------------------------
ecbuild_add_test( TARGET  test_soca_geometry
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_state
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_modelaux
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_model
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_increment
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_errorcovariance
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/static_error_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_linearmodel
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linearmodel_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_balance
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/balance.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerrfilt
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/bkgerrfilt.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerrsoca
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/bkgerrsoca.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerrgodas
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/bkgerrgodas.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_vertconv
                  MPI ${MPI_SOCA_TESTS}
                  SOURCES executables/TestVariableChange.cc
                  ARGS    "testinput/vertconv.yml"
                  LIBS    soca )


#-------------------------------------------------------------------------------
# Tests using the main SOCA applications
#  TEST output stream is compared against existing testoutput/ files
#-------------------------------------------------------------------------------
soca_add_test( NAME ensvariance
               EXE  soca_ensvariance.x)

soca_add_test( NAME parameters_bump_loc
               EXE  soca_parameters.x
               NOCOMPARE)

#soca_add_test( NAME parameters_bump_cov_lct
#               EXE  soca_parameters.x)

#soca_add_test( NAME parameters_bump_cov_nicas
#               EXE  soca_parameters.x )

# Dirac tests
#soca_add_test( NAME dirac_bump_cov
#               EXE  soca_dirac.x )

soca_add_test( NAME dirac_soca_cov
               EXE  soca_dirac.x )

soca_add_test( NAME dirac_socahyb_cov
               EXE  soca_dirac.x )

#ecbuild_add_test( TARGET  test_soca_localization
#                  BOOST
#                  SOURCES executables/TestLocalization.cc
#                  ARGS    "testinput/localization_test.yml"
#                  MPI     1
#                  LIBS    soca )

#ecbuild_add_test( TARGET test_soca_obserrorcov
#                  BOOST
#                  SOURCES executables/TestObsErrorCovariance.cc
#                  ARGS "testinput/obserrorcovariance_test.yml"
#                  LIBS    soca )
#

# h(x) executables
soca_add_test( NAME hofx3d
               EXE  soca_hofx3d.x)

soca_add_test( NAME hofx
               EXE  soca_hofx.x)

soca_add_test( NAME enshofx
               EXE  soca_enshofx.x)

# variational methods
soca_add_test( NAME 3dvar_soca
               EXE  soca_3dvar.x)

soca_add_test( NAME 3dvar_godas
               EXE  soca_3dvar.x)

soca_add_test( NAME 3dvarfgat
               EXE  soca_3dvar.x)

soca_add_test( NAME 3dhyb
               EXE  soca_3dvar.x)

soca_add_test( NAME 3dhybfgat
               EXE  soca_3dvar.x)

#soca_add_test( NAME 4dhyb
#               EXE  soca_3dvar.x)

# restart checkpointing
soca_add_test( NAME checkpointmodel
               EXE  soca_checkpoint_model.x)
