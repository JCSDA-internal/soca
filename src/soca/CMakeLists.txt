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

target_compile_features( soca PUBLIC cxx_std_11 )

target_link_libraries( soca PUBLIC NetCDF::NetCDF_Fortran )
target_link_libraries( soca PUBLIC fckit )
target_link_libraries( soca PUBLIC atlas )
target_link_libraries( soca PUBLIC atlas_f )
target_link_libraries( soca PUBLIC fms_r8 )
target_link_libraries( soca PUBLIC gsw )
target_link_libraries( soca PUBLIC mom6 )
target_link_libraries( soca PUBLIC oops )
target_link_libraries( soca PUBLIC saber )
target_link_libraries( soca PUBLIC ioda )
target_link_libraries( soca PUBLIC ufo )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${DEPEND_LIB_ROOT}/mod )
# Add source code in the subdirectories
add_subdirectory(AnalyticInit)
add_subdirectory(Covariance)
add_subdirectory(Fields)
add_subdirectory(Geometry)
add_subdirectory(GeometryIterator)
add_subdirectory(Increment)
add_subdirectory(Interpolator)
add_subdirectory(LinearVariableChange)
add_subdirectory(Model)
add_subdirectory(ModelBias)
add_subdirectory(ObsLocalization)
add_subdirectory(State)
add_subdirectory(Utils)
add_subdirectory(VariableChange)

separate_arguments(ESMF_F90LINKPATHS NATIVE_COMMAND ${ESMF_F90LINKPATHS})
set(tmp "")
foreach (ITEM ${ESMF_F90LINKPATHS})
   string(REGEX REPLACE "-L" "" ITEM "${ITEM}")
   list(APPEND tmp ${ITEM})
endforeach()
set(ESMF_F90LINKPATHS ${tmp})
target_link_libraries( soca PUBLIC ufs)
target_link_libraries( soca PUBLIC NetCDF::NetCDF_Fortran)
target_link_libraries( soca PUBLIC nemsio::nemsio)
target_link_libraries( soca PUBLIC w3emc::w3emc_d)
target_link_libraries( soca PUBLIC w3nco::w3nco_d)
target_link_libraries( soca PUBLIC sp::sp_d)
target_link_libraries( soca PUBLIC bacio::bacio_4)
target_link_libraries( soca PUBLIC FMS::fms_r8)
target_link_libraries( soca PUBLIC esmf)
target_link_directories( soca PUBLIC ${ESMF_F90LINKPATHS} )
target_compile_definitions( soca PUBLIC "-DFRONT_MOM6=mom_cap_mod")
target_compile_definitions(soca PUBLIC UFS_BUILD)
add_dependencies(soca FMS )
add_dependencies(soca ufs-weather-model )

## Fortran modules
set(MODULE_DIR ${PROJECT_NAME}/module)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR}>)
