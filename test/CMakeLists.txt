list( APPEND mom5cice5_test_input
  testinput/dirac.nicas.json
  testinput/interfaces.json
  testinput/test.json
  compare.sh
)

list( APPEND mom5cice5_testoutput
  testoutput/dirac.nicas.test
)

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${mom5cice5_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   mom5cice5_test_scripts
                       SOURCES_PACK
                       ${mom5cice5_test_input}
                     )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

#####################################################################

ecbuild_add_test( TARGET  test_mom5cice5_geometry
                  BOOST
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    mom5cice5 )

ecbuild_add_test( TARGET  test_mom5cice5_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    mom5cice5 )

ecbuild_add_test( TARGET  test_mom5cice5_increment
                  BOOST
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/interfaces.json"
                  LIBS    mom5cice5 )

ecbuild_add_test( TARGET test_mom5cice5_dirac_nicas
                  TYPE SCRIPT
                  COMMAND "compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/mom5cice5_dirac.x testinput/dirac.nicas.json"
                       testoutput/dirac.nicas.test
                  DEPENDS mom5cice5_dirac.x )

