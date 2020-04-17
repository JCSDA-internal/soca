# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# only the source code in this directory, all other files
# are handled through the subdirectory's respective CMakeLists.txt
list( APPEND soca_src_files
  Traits.h
  Fortran.h
)

# a convenience function to make adding source files from subdirectories easier
# (to keep from having to use the full path every time)
function(soca_target_sources)
  foreach( src IN LISTS ARGN)
    target_sources(soca PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${src})
  endforeach()
endfunction()

# The main soca library... behold its glory
ecbuild_add_library( TARGET   soca
                     SOURCES ${soca_src_files}
                     LIBS     ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES}
                              eckit eckit_mpi fckit oops saber ioda ufo fms mom6
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${SOCA_LINKER_LANGUAGE}
                    )

# Add source code in the subdirectories
add_subdirectory(Covariance)
add_subdirectory(Fields)
add_subdirectory(Geometry)
add_subdirectory(GeometryIterator)
add_subdirectory(GetValues)
add_subdirectory(Increment)
add_subdirectory(LinearModel)
add_subdirectory(Localization)
add_subdirectory(Model)
add_subdirectory(ModelBias)
add_subdirectory(State)
add_subdirectory(Transforms)
add_subdirectory(Utils)
