list( APPEND soca_test_input
  testinput/forecast.yml 
  testinput/state_test.yml
  testinput/model_test.yml
  testinput/linearmodel_test.yml
  testinput/geometry_test.yml
  testinput/increment_test.yml
  testinput/hofx-adt.yml
  testinput/3dvar_test.yml
  testinput/3dvarfgat_test.yml  
  testinput/3dhybridbump.yml  
  testinput/balance.yml
  testinput/bkgerr.yml
  testinput/vertconv.yml
  testinput/genenspert.yml
  testinput/static_SocaError_init.yml
  testinput/static_SocaError_test.yml
  compare.sh
  run.sh
)

list( APPEND soca_testoutput
  testoutput/forecast.test
  testoutput/hofx.test  
  testoutput/3dvar.test  
)

list( APPEND soca_model_param
  Data/360x210x63/input.nml
  Data/360x210x63/MOM_input
  Data/360x210x63/MOM_override
  Data/360x210x63/Vertical_coordinate.nc
  Data/360x210x63/ice.bkgerror.nc
  Data/360x210x63/ocn.bkgerror.nc
  Data/rossrad.dat
)

list( APPEND soca_model_restarts
  Data/360x210x63/INPUT/MOM.res.nc
  Data/360x210x63/INPUT/ice_model.res.nc
  Data/360x210x63/INPUT/ocean_hgrid.nc
  Data/360x210x63/INPUT/topog.nc
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
    message(${FILENAME})
    message(${CMAKE_CURRENT_SOURCE_DIR})    
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

#Link to MOM's resource files
foreach(FILENAME ${soca_model_param})
  get_filename_component(filename ${FILENAME} NAME )
  message("filename="${filename})
     message("FILENAME="${FILENAME})     
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

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
     message("filename="${filename})
     message("FILENAME="${FILENAME})
     message(${CMAKE_CURRENT_BINARY_DIR})
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )	   
endforeach(FILENAME)

# Create directory for BUMP output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)


#####################################################################

# Tests that create data other tests will use

# Test of read/write restart, model advance = Id 
ecbuild_add_test( TARGET test_soca_forecast
                  TYPE SCRIPT
		  #MPI     4
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 2 ${CMAKE_BINARY_DIR}/bin/soca_forecast.x testinput/forecast.yml"
                  DEPENDS soca_forecast.x )

# Initialize static error covariance model defined in soca, used for 3DVAR
ecbuild_add_test( TARGET test_soca_socaerror_init
                  TYPE SCRIPT
		  #MPI     4
		  COMMAND "./run.sh"
                  ARGS "mpirun -np 2 ${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x testinput/static_SocaError_init.yml"
                  DEPENDS soca_staticbinit.x )

#ecbuild_add_test( TARGET test_soca_socaerror_serial_init
#                  TYPE SCRIPT
#		  #MPI     4
#		  COMMAND "./run.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x testinput/static_SocaError_init.yml"
#                  DEPENDS soca_staticbinit.x )


# Initialize static error covariance model for the generation of perturbation
#ecbuild_add_test( TARGET test_soca_enspert_B_init
#                  TYPE SCRIPT
#		  COMMAND "./run.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x testinput/enspert_B_init.yml"
#                  DEPENDS soca_staticbinit.x )

ecbuild_add_test( TARGET test_soca_enspert
                  TYPE SCRIPT
		  #MPI     4
		  COMMAND "./run.sh"
                  ARGS "mpirun -np 2 ${CMAKE_BINARY_DIR}/bin/soca_enspert.x testinput/genenspert.yml"
                  DEPENDS soca_enspert.x )
		
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
		  #MPI     4
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_linearmodel
                  BOOST
		  #MPI     4		  
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linearmodel_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_state
                  BOOST
		  #MPI     4
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
		  #MPI     4
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_increment
                  BOOST
		  #MPI     4 		  
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_errorcovariance
                  BOOST
		  #MPI     4
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/static_SocaError_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_balance
                  BOOST
		  #MPI     4
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/balance.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerr
                  BOOST
		  #MPI     4
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/bkgerr.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_vertconv
                  BOOST
		  #MPI     4
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/vertconv.yml"
                  LIBS    soca )		
		
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
                  TYPE SCRIPT
		  #MPI     4
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 2 ${CMAKE_BINARY_DIR}/bin/soca_hofx.x testinput/hofx-adt.yml"             
                  DEPENDS soca_hofx.x )

ecbuild_add_test( TARGET test_soca_3dvar
                  TYPE SCRIPT
		  #MPI     4
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 2 ${CMAKE_BINARY_DIR}/bin/soca_3dvar.x testinput/3dvar_test.yml"
		         testoutput/3dvar.test
                  DEPENDS soca_3dvar.x )

ecbuild_add_test( TARGET test_soca_3dvarfgat
                  TYPE SCRIPT
		  #MPI     4
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 2 ${CMAKE_BINARY_DIR}/bin/soca_3dvar.x testinput/3dvarfgat_test.yml"
		         testoutput/3dvarfgat.test
                  DEPENDS soca_3dvar.x )
		
ecbuild_add_test( TARGET test_soca_3dhybrid
                  TYPE SCRIPT
		  #MPI     4
                  COMMAND "./compare.sh"
                  ARGS "$mpirun -np 2 {CMAKE_BINARY_DIR}/bin/soca_3dvar.x testinput/3dhybridbump.yml"
                  DEPENDS soca_3dvar.x )
		  
#ecbuild_add_test( TARGET test_soca_dirac
#                  TYPE SCRIPT
#                  MPI     1		  
#                  COMMAND "./compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_dirac.x testinput/nicas_test.yml"
#                       testoutput/dirac.test
#                  DEPENDS soca_3dvar.x )		  
  
#ecbuild_add_test( TARGET test_soca_parameters_bump_cov
#                  TYPE SCRIPT
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_parameters.x testinput/parametersbump.cov.yml" )
