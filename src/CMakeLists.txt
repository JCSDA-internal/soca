# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( soca )
add_subdirectory( mains )

# find cpplint.py, if not found assume it will be built by oops in this bundle
find_file(CPPLINT_PY cpplint.py)
if(NOT CPPLINT_PY)
  set(CPPLINT_PY ${CMAKE_BINARY_DIR}/bin/cpplint.py)
endif()

ecbuild_add_test( TARGET soca_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CPPLINT_PY}
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
