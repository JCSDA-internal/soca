list( APPEND soca_test_input
  testinput/interfaces.json  
  testinput/truth.json
  testoutput/forecast.test  
  testinput/makeobs3d.yml
  testinput/makeobs3d.json
  testinput/state_test.json    
  testinput/geometry_test.json
  testinput/increment_test.json
  compare.sh
)

list( APPEND soca_testoutput
  #testoutput/dirac.nicas.test
  testoutput/forecast.test
  testoutput/truth.test
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

ecbuild_add_resources( TARGET   soca_test_scripts
                       SOURCES_PACK
                       ${soca_test_input}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

ecbuild_add_test( TARGET test_soca_truth
                  TYPE SCRIPT
		  MPI     2  
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_forecast.x testinput/truth.json"
                       testoutput/truth.test
                  DEPENDS soca_forecast.x )

ecbuild_add_test( TARGET test_soca_forecast
                  TYPE SCRIPT
		  MPI     2		  		  
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_forecast.x testinput/forecast.yml"
                       testoutput/forecast.test
                  DEPENDS soca_forecast.x )

#ecbuild_add_test( TARGET test_soca_makeobs3d
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_makeobs.x testinput/makeobs3d.json"
#                       testoutput/makeobs3d.test
#                  DEPENDS soca_makeobs.x )
		  
# Tests of class interfaces

ecbuild_add_test( TARGET  test_soca_geometry
                  BOOST
		  MPI     2
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_test.json"
                  LIBS    soca )

#ecbuild_add_test( TARGET  test_soca_locations
#                  BOOST
#                  SOURCES executables/TestLocations.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )

#ecbuild_add_test( TARGET  test_soca_variables
#                  BOOST
#		  MPI     2		  
#                  SOURCES executables/TestVariables.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_state
                  BOOST
		  MPI     2
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state_test.json"
                  LIBS    soca )
		  
#ecbuild_add_test( TARGET  test_soca_modelaux
#                  BOOST
#                  SOURCES executables/TestModelAuxControl.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )
		  
#ecbuild_add_test( TARGET  test_soca_model
#                  BOOST
#                  SOURCES executables/TestModel.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_increment
                  BOOST
		  MPI     2 		  
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_test.json"
                  LIBS    soca )

#ecbuild_add_test( TARGET  test_soca_errorcovariance
#                  BOOST
#                  SOURCES executables/TestErrorCovariance.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )

#ecbuild_add_test( TARGET  test_soca_obsvector
#                  BOOST
#                  SOURCES executables/TestObsVector.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )

#ecbuild_add_test( TARGET  test_soca_obsoperator
#                  BOOST
#                  SOURCES executables/TestObsOperator.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )

#ecbuild_add_test( TARGET  test_soca_linobsoperator
#                  BOOST
#                  SOURCES executables/TestLinearObsOperator.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )

#ecbuild_add_test( TARGET test_soca_observationspace
#                  BOOST
#                  SOURCES executables/TestObservationSpace.cc
#                  ARGS "testinput/interfaces.json"
#                  LIBS    soca )

#ecbuild_add_test( TARGET  test_soca_obsaux
#                  BOOST
#                  SOURCES executables/TestObsAuxControl.cc
#                  ARGS    "testinput/interfaces.json"
#                  LIBS    soca )
		  
#ecbuild_add_test( TARGET test_soca_dirac_nicas
#                  TYPE SCRIPT
#                  COMMAND "compare.sh"
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_dirac.x testinput/dirac.nicas.json"
#                       testoutput/dirac.nicas.test
#                  DEPENDS soca_dirac.x )

