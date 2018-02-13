list( APPEND soca_test_input
  testinput/truth.yml
  testinput/forecast.yml 
  testinput/makeobs3d.yml
  testinput/state_test.yml
  testinput/hofx_test.yml
  testinput/model_test.yml
  testinput/linearmodel_test.yml
  testinput/geometry_test.yml
  testinput/increment_test.yml
  testinput/errorcovariance_test.yml
  testinput/obserrorcovariance_test.yml  
  testinput/3dvar_test.yml
  compare.sh
)

list( APPEND soca_testoutput
  #testoutput/dirac.nicas.test
  testoutput/forecast.test
  testoutput/truth.test
  testoutput/hofx.test  
  testoutput/makeobs3d.test
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

# Link to MOM6-SIS2 scratch dir files
file( GLOB mom6_scratch_files ${MOM6_SCRATCH_DIR}/* )
foreach(FILENAME ${mom6_scratch_files})
     get_filename_component(filename ${FILENAME} NAME )
     message("Symbolic link of "${filename}" to:"${FILENAME})     
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${MOM6_SCRATCH_DIR}/*
           ${CMAKE_CURRENT_BINARY_DIR}/ )

ecbuild_add_resources( TARGET   soca_test_scripts
                       SOURCES_PACK
                       ${soca_test_input}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${UFO_DATA_DIR}/seaice_obs.nc
           ${CMAKE_CURRENT_BINARY_DIR}/Data/seaice_obs.nc )


#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

ecbuild_add_test( TARGET test_soca_truth
                  TYPE SCRIPT
		  MPI     2  
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_forecast.x testinput/truth.yml"
                       testoutput/truth.test
                  DEPENDS soca_forecast.x )

ecbuild_add_test( TARGET test_soca_forecast
                  TYPE SCRIPT
		  MPI     2		  		  
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_forecast.x testinput/forecast.yml"
                       testoutput/forecast.test
                  DEPENDS soca_forecast.x )

ecbuild_add_test( TARGET test_soca_makeobs3d
                  TYPE SCRIPT
		  MPI     1  		  		  
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_makeobs.x testinput/makeobs3d.yml"
                       testoutput/makeobs3d.test
                  DEPENDS soca_makeobs.x )
		  
# Tests of class interfaces

ecbuild_add_test( TARGET  test_soca_geometry
                  BOOST
		  MPI     2
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_linearmodel
                  BOOST
		  MPI     2		  
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
		  MPI     2
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_increment
                  BOOST
		  MPI     2 		  
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/errorcovariance_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_obserrorcov
                  BOOST
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/obserrorcovariance_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_hofx
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_hofx.x testinput/hofx_test.yml"
                       testoutput/hofx.test
                  DEPENDS soca_hofx.x )

ecbuild_add_test( TARGET test_soca_3dvar
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_3dvar.x testinput/3dvar_test.yml"
                       testoutput/3dvar.test
                  DEPENDS soca_3dvar.x )
		  
ecbuild_add_test( TARGET test_soca_dirac
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_dirac.x testinput/dirac_test.json"
                       testoutput/3dvar.test
                  DEPENDS soca_3dvar.x )		  
