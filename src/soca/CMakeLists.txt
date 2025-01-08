# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

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
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${SOCA_LINKER_LANGUAGE}
                    )

target_include_directories(soca PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
  $<INSTALL_INTERFACE:include/soca>)

target_link_libraries( soca PUBLIC NetCDF::NetCDF_Fortran )
target_link_libraries( soca PUBLIC NetCDF::NetCDF_C )
target_link_libraries( soca PUBLIC fckit )
target_link_libraries( soca PUBLIC atlas )
target_link_libraries( soca PUBLIC atlas_f )
target_link_libraries( soca PUBLIC fms )
target_link_libraries( soca PUBLIC gsw )
target_link_libraries( soca PUBLIC mom6 )
target_link_libraries( soca PUBLIC icepack_limited ) # custom soca build with only some functions built
target_link_libraries( soca PUBLIC oops )
target_link_libraries( soca PUBLIC saber )
target_link_libraries( soca PUBLIC ioda )
target_link_libraries( soca PUBLIC ufo )
target_link_libraries( soca PUBLIC vader )


if(Torch_FOUND)
  target_link_libraries( soca PUBLIC
    "${Torch_ROOT}/lib/libtorch.so"          # For the Torch core library
    "${Torch_ROOT}/lib/libtorch_cpu.so"      # For the CPU backend
    "${Torch_ROOT}/lib/libc10.so"            # Core library
    )
endif()

# Add source code in the subdirectories
add_subdirectory(AnalyticInit)
add_subdirectory(Covariance)
add_subdirectory(Fields)
add_subdirectory(Geometry)
add_subdirectory(GeometryIterator)
add_subdirectory(Increment)
add_subdirectory(LinearVariableChange)
if(Torch_FOUND)
  add_subdirectory(MLBalance)
endif()
add_subdirectory(ModelBias)
add_subdirectory(ModelData)
add_subdirectory(ObsLocalization)
add_subdirectory(SaberBlocks)
add_subdirectory(State)
add_subdirectory(Utils)
add_subdirectory(VariableChange)

## Fortran modules
set(MODULE_DIR ${PROJECT_NAME}/module)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR}>)
