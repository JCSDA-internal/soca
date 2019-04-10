list( APPEND soca_test_input
  testinput/checkpointmodel.yml  
  testinput/ensvariance.yml
  testinput/forecast.yml
  testinput/forecast_mom6.yml
  testinput/state_test.yml
  testinput/model_test.yml
  testinput/linearmodel_test.yml
  testinput/geometry_test.yml
  testinput/increment_test.yml
  testinput/hofx.yml
  testinput/3dvar_soca.yml
  testinput/3dvar_godas.yml  
  testinput/3dvarfgat_mom6.yml  
  testinput/4dhyb.yml  
  testinput/balance.yml
  testinput/bkgerr.yml
  testinput/bkgerrfilt.yml
  testinput/bkgerrsoca.yml
  testinput/bkgerrgodas.yml
  testinput/vertconv.yml
  testinput/genenspert.yml
  testinput/ensvariance.yml
  testinput/static_SocaError_init.yml
  testinput/static_SocaLoc_init.yml
  testinput/static_SocaError_test.yml
  testinput/parameters_bump_cov_lct.yml
  testinput/parameters_bump_cov_nicas.yml
  testinput/dirac_bump_cov.yml
  testinput/dirac_soca_cov.yml
  compare.sh
)

list( APPEND soca_testoutput
  testoutput/forecast.test
  testoutput/hofx.test
  testoutput/3dvar.test
  testoutput/genenspert.test
  testoutput/parameters_bump_cov_lct.test
  testoutput/parameters_bump_cov_nicas.test
  testoutput/dirac_bump_cov.test
  testoutput/dirac_soca_cov.test
)

list( APPEND soca_model_param
  Data/360x210x75/input.nml
  Data/360x210x75/MOM_input
  Data/360x210x75/MOM_override
  Data/360x210x75/data_table
  Data/360x210x75/diag_table
  Data/360x210x75/field_table
  Data/360x210x75/ice.bkgerror.nc	
  Data/360x210x75/ocn.bkgerror.nc
  Data/rossrad.dat
  Data/godas_sst_bgerr.nc  
)

list( APPEND soca_model_restarts
  Data/360x210x75/INPUT/MOM.res.nc
  Data/360x210x75/INPUT/ice_model.res.nc
  Data/360x210x75/INPUT/layer_coord.nc
  Data/360x210x75/INPUT/hycom1_75_800m.nc
  Data/360x210x75/INPUT/ocean_hgrid.nc
  Data/360x210x75/INPUT/topog.nc
)


list( APPEND soca_model_forcing
  Data/360x210x75/INPUT/geothermal_heating_cm2g.nc
  Data/360x210x75/INPUT/seawifs_1998-2006_GOLD_smoothed_2X.nc
  Data/360x210x75/INPUT/tideamp.nc
  Data/360x210x75/INPUT/topog.nc
  Data/360x210x75/INPUT/gustiness_qscat.nc
  Data/360x210x75/INPUT/ocean_forcing_daily.nc
  Data/360x210x75/INPUT/ocean_precip_monthly.nc
  Data/360x210x75/INPUT/sgs_h2.nc
)

list( APPEND soca_obs
  Data/Obs/Jason-2-2018-04-15-soca.nc
)

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${soca_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${soca_testoutput})
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

# MOM's forcing files
foreach(FILENAME ${soca_model_forcing})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/INPUT/${filename} )
endforeach(FILENAME)

# MOM's output directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RESTART)

ecbuild_add_resources( TARGET   soca_test_scripts
                       SOURCES_PACK
                       ${soca_test_input}
                       )
		     
# Link to IODA Marine obs
file(GLOB ioda_obs "${ioda_SOURCE_DIR}/test/testinput/marine/*.nc*")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
foreach(FILENAME ${ioda_obs})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )	   
endforeach(FILENAME)

# Link to SOCA Marine obs
foreach(FILENAME ${soca_obs})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )	   
endforeach(FILENAME)

# Create directory for BUMP output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)


#####################################################################
# NOTE: MPI defaults can be overriden by specifying -DMPIEXEC_MAX_NUMPROCS=?
#   when running ecbuild


# Tests that create data other tests will use

# Test of read/write restart, model advance = Id
ecbuild_add_test( TARGET test_soca_forecast
                  SOURCES ../mains/Forecast.cc
                  LIBS soca
                  ARGS "testinput/forecast.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS})

# Test mom6 forecast
ecbuild_add_test( TARGET test_soca_forecast_mom6
                  SOURCES ../mains/Forecast.cc
                  LIBS soca
                  ARGS "testinput/forecast_mom6.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

# Initialize static error covariance model defined in soca, used for 3DVAR
ecbuild_add_test( TARGET test_soca_socaerror_init
                  SOURCES ../mains/StaticBInit.cc
                  LIBS soca
                  ARGS "testinput/static_SocaError_init.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_socaloc_init
                  SOURCES ../mains/StaticBInit.cc
                  LIBS soca
                  ARGS "testinput/static_SocaLoc_init.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_socaerror_serial_init
                  SOURCES ../mains/StaticBInit.cc
                  LIBS soca
                  ARGS "testinput/static_SocaError_init.yml"
                  MPI 1 )

ecbuild_add_test( TARGET test_soca_enspert
                  TYPE SCRIPT
                  COMMAND "./compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_enspert.x testinput/genenspert.yml"
                       testoutput/genenspert.test
                  DEPENDS soca_enspert.x
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_ensvariance
                  SOURCES ../mains/EnsVariance.cc
                  LIBS soca
                  ARGS "testinput/ensvariance.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

# Test below can take a few hours to complete. Commented out for now
#ecbuild_add_test( TARGET test_soca_covar_lengthscale
#                  TYPE SCRIPT
#		  COMMAND "./run.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x testinput/static_scale_binit.yml"
#                  DEPENDS soca_staticbinit.x )

#ecbuild_add_test( TARGET test_soca_binit
#                  TYPE SCRIPT
#		  COMMAND "./run.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x testinput/static_binit.yml"
#                  DEPENDS soca_staticbinit.x )


# Tests of class interfaces

ecbuild_add_test( TARGET  test_soca_geometry
                  BOOST
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_linearmodel
                  BOOST
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linearmodel_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_modelaux
                  BOOST
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_model
                  BOOST
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/static_SocaError_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_balance
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/balance.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerrfilt
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/bkgerrfilt.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerrsoca
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/bkgerrsoca.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerrgodas
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/bkgerrgodas.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_vertconv
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/vertconv.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_dirac_soca_cov
                  TYPE SCRIPT
                  COMMAND "./compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_dirac.x testinput/dirac_soca_cov.yml"
                       testoutput/dirac_soca_cov.test
                  DEPENDS soca_dirac.x 
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_parameters_bump_cov_lct
                  TYPE SCRIPT
                  COMMAND "./compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_parameters.x testinput/parameters_bump_cov_lct.yml"
                       testoutput/parameters_bump_cov_lct.test
                  DEPENDS soca_parameters.x 
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_parameters_bump_cov_nicas
                  TYPE SCRIPT
                  COMMAND "./compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_parameters.x testinput/parameters_bump_cov_nicas.yml"
                       testoutput/parameters_bump_cov_nicas.test
                  DEPENDS soca_parameters.x 
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_dirac_bump_cov
                  TYPE SCRIPT
                  COMMAND "./compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_dirac.x testinput/dirac_bump_cov.yml"
                       testoutput/dirac_bump_cov.test
                  DEPENDS soca_dirac.x 
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

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

ecbuild_add_test( TARGET test_soca_hofx
                  SOURCES ../mains/HofX.cc
                  LIBS soca
                  ARGS "testinput/hofx.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_3dvarsoca
                  SOURCES ../mains/3DVar.cc
                  LIBS soca
                  ARGS "testinput/3dvar_soca.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_checkpointmodel
                  SOURCES ../mains/CheckpointModel.cc
                  LIBS soca
                  ARGS "testinput/checkpointmodel.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_3dvargodas
                  SOURCES ../mains/3DVar.cc
                  LIBS soca
                  ARGS "testinput/3dvar_godas.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )		

ecbuild_add_test( TARGET test_soca_3dvarfgat_mom6
                  SOURCES ../mains/3DVar.cc
                  LIBS soca
                  ARGS "testinput/3dvarfgat_mom6.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )

ecbuild_add_test( TARGET test_soca_4dhybrid
                  SOURCES ../mains/3DVar.cc
                  LIBS soca
                  ARGS "testinput/4dhyb.yml"
                  MPI ${MPIEXEC_MAX_NUMPROCS} )
