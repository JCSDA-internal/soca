# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( soca )
add_subdirectory( mains )

# find cpplint.py provided by oops
find_file(CPPLINT_PY oops_cpplint.py HINTS ${oops_BINDIR})

ecbuild_add_test( TARGET soca_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CPPLINT_PY}
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
