list( APPEND soca_test_input
  compare.sh
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
  testinput/ensvariance.yml
  testinput/forecast_identity.yml
  testinput/forecast_identity_cice.yml
  testinput/forecast_mom6.yml
  testinput/genenspert.yml
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
  testinput/static_SocaError_init.yml
  testinput/static_SocaError_test.yml
  testinput/vertconv.yml
)

list( APPEND soca_test_output
  testoutput/3dhyb.test
  testoutput/3dhybfgat.test
  testoutput/3dvar_godas.test
  testoutput/3dvar_soca.test
  testoutput/3dvarfgat.test
  testoutput/checkpointmodel.test
  testoutput/dirac_bump_cov.test
  testoutput/dirac_soca_cov.test
  testoutput/dirac_socahyb_cov.test
  testoutput/enshofx.test
  testoutput/ensvariance.test
  testoutput/forecast_identity.test
  testoutput/forecast_identity_cice.test
  testoutput/forecast_mom6.test
  testoutput/genenspert.test
  testoutput/gridgen.test
  testoutput/hofx.test
  testoutput/hofx3d.test
  testoutput/parameters_bump_cov_nicas.test
  testoutput/parameters_bump_loc.test
)

list( APPEND soca_model_param
  Data/90x43x25/diag_table
  Data/90x43x25/field_table
  Data/90x43x25/ice.bkgerror.nc
  Data/90x43x25/input.nml
  Data/90x43x25/MOM_input
  Data/90x43x25/ocn.bkgerror.nc
  Data/godas_sst_bgerr.nc
  Data/rossrad.dat
)

list( APPEND soca_model_restarts
  Data/90x43x25/INPUT/forcing_daily.nc
  Data/90x43x25/INPUT/forcing_monthly.nc
  Data/90x43x25/INPUT/hycom1_25.nc
  Data/90x43x25/INPUT/ice_model.res.nc
  Data/90x43x25/INPUT/cice.res.nc
  Data/90x43x25/INPUT/layer_coord25.nc
  Data/90x43x25/INPUT/MOM.res.nc
  Data/90x43x25/INPUT/sfc.res.nc
  Data/90x43x25/INPUT/ocean_hgrid.nc
  Data/90x43x25/INPUT/ocean_topog.nc
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
foreach(FILENAME ${soca_test_input} ${soca_test_output})
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
set(MPI_SOCA_TESTS 2)

#-------------------------------------------------------------------------------
# Tests that create data other tests will use
#-------------------------------------------------------------------------------

# Test of grid generation
ecbuild_add_test( TARGET test_soca_generate_grid
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_gridgen.x
                      testinput/gridgen.yml"
                      testoutput/gridgen.test
                  DEPENDS soca_gridgen.x)

# Test of read/write restart for MOM6&SIS2, model advance = Id
ecbuild_add_test( TARGET test_soca_forecast_identity
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_forecast.x
                      testinput/forecast_identity.yml"
                      testoutput/forecast_identity.test
                  DEPENDS soca_forecast.x)

# Test of read/write restart for MOM6&CICE, model advance = Id
ecbuild_add_test( TARGET test_soca_forecast_identity_cice
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_forecast.x
                      testinput/forecast_identity_cice.yml"
                      testoutput/forecast_identity_cice.test
                  DEPENDS soca_forecast.x)

# Test mom6 forecast
ecbuild_add_test( TARGET test_soca_forecast_mom6
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_forecast.x
                      testinput/forecast_mom6.yml"
                      testoutput/forecast_mom6.test
                  DEPENDS soca_forecast.x)

# Initialize static error covariance model defined in soca, used for 3DVAR
ecbuild_add_test( TARGET test_soca_socaerror_init
                  MPI ${MPI_SOCA_TESTS}
                  TYPE EXE
                  COMMAND ${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x
                  ARGS testinput/static_SocaError_init.yml
                  DEPENDS soca_staticbinit.x)

ecbuild_add_test( TARGET test_soca_enspert
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_enspert.x
                      testinput/genenspert.yml"
                      testoutput/genenspert.test
                  DEPENDS soca_enspert.x)

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
                  ARGS    "testinput/static_SocaError_test.yml"
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

ecbuild_add_test( TARGET test_soca_ensvariance
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_ensvariance.x
                      testinput/ensvariance.yml"
                      testoutput/ensvariance.test
                  DEPENDS soca_ensvariance.x)

# The following two tests are slow, how can we speed it up?
# Use another model grid configuration with even fewer grid cells?
ecbuild_add_test( TARGET test_soca_parameters_bump_loc
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_parameters.x
                      testinput/parameters_bump_loc.yml"
                      testoutput/parameters_bump_loc.test
                  DEPENDS soca_parameters.x)

ecbuild_add_test( TARGET test_soca_parameters_bump_cov_lct
                  ENABLED OFF
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_parameters.x
                      testinput/parameters_bump_cov_lct.yml"
                      testoutput/parameters_bump_cov_lct.test
                  DEPENDS soca_parameters.x)

ecbuild_add_test( TARGET test_soca_parameters_bump_cov_nicas
                  ENABLED OFF
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_parameters.x
                      testinput/parameters_bump_cov_nicas.yml"
                      testoutput/parameters_bump_cov_nicas.test
                  DEPENDS soca_parameters.x)

ecbuild_add_test( TARGET test_soca_dirac_bump_cov
                  ENABLED OFF
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_dirac.x
                      testinput/dirac_bump_cov.yml"
                      testoutput/dirac_bump_cov.test
                  DEPENDS soca_dirac.x )

ecbuild_add_test( TARGET test_soca_dirac_soca_cov
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_dirac.x
                      testinput/dirac_soca_cov.yml"
                      testoutput/dirac_soca_cov.test
                  DEPENDS soca_dirac.x )

ecbuild_add_test( TARGET test_soca_dirac_socahyb_cov
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_dirac.x
                      testinput/dirac_socahyb_cov.yml"
                      testoutput/dirac_socahyb_cov.test
                  DEPENDS soca_dirac.x )

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
ecbuild_add_test( TARGET test_soca_hofx3d
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_hofx3d.x
                      testinput/hofx3d.yml"
                      testoutput/hofx3d.test
                  DEPENDS soca_hofx3d.x)

ecbuild_add_test( TARGET test_soca_hofx
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_hofx.x
                      testinput/hofx.yml"
                      testoutput/hofx.test
                  DEPENDS soca_hofx.x)

ecbuild_add_test( TARGET test_soca_enshofx
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_enshofx.x
                      testinput/enshofx.yml"
                      testoutput/enshofx.test
                  DEPENDS soca_hofx.x)

ecbuild_add_test( TARGET test_soca_3dvarsoca
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_3dvar.x
                      testinput/3dvar_soca.yml"
                      testoutput/3dvar_soca.test
                  DEPENDS soca_3dvar.x)

ecbuild_add_test( TARGET test_soca_3dvargodas
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_3dvar.x
                      testinput/3dvar_godas.yml"
                      testoutput/3dvar_godas.test
                  DEPENDS soca_3dvar.x)

ecbuild_add_test( TARGET test_soca_3dhyb
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_3dvar.x
                      testinput/3dhyb.yml"
                      testoutput/3dhyb.test
                  DEPENDS soca_3dvar.x)

ecbuild_add_test( TARGET test_soca_3dhybfgat
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_3dvar.x
                      testinput/3dhybfgat.yml"
                      testoutput/3dhybfgat.test
                  DEPENDS soca_3dvar.x)

ecbuild_add_test( TARGET test_soca_checkpointmodel
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_checkpoint_model.x
                      testinput/checkpointmodel.yml"
                      testoutput/checkpointmodel.test
                  DEPENDS soca_checkpoint_model.x)

ecbuild_add_test( TARGET test_soca_3dvarfgat
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_3dvar.x
                      testinput/3dvarfgat.yml"
                      testoutput/3dvarfgat.test
                  DEPENDS soca_3dvar.x)

ecbuild_add_test( TARGET test_soca_4dhybrid
                  ENABLED OFF
                  MPI ${MPI_SOCA_TESTS}
                  TYPE SCRIPT
                  COMMAND ./compare.sh
                  ARGS
                      "${CMAKE_BINARY_DIR}/bin/soca_3dvar.x
                      testinput/4dhyb.yml"
                      testoutput/4dhyb.test
                  DEPENDS soca_3dvar.x)
