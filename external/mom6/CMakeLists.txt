# This file was modified from the original wrapper in ufs-weather-model
# https://github.com/ufs-community/ufs-weather-model/tree/develop/MOM6-interface

### MOM6 Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Waliasing -fcray-pointer -fconvert=big-endian -ffree-line-length-none -fno-range-check -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -sox")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  set(CMAKE_Fortran_LINK_FLAGS "")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MOM6_Extra_FORTRAN_FLAGS}")

# User option to build MOM6 solo executable
set(MOM6SOLO ON CACHE BOOL "Build MOM6 solo")

# include the rats nest of mom6 files
include("mom6_files.cmake")

### Use common object library for building target library
add_library(mom6 STATIC
    ${mom6_src_files}
    ${mom6_stochastic_stub_files} ) ## TODO include stochastic package at some point?
set_target_properties(mom6 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(mom6 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
    $<INSTALL_INTERFACE:mod>)
target_include_directories(mom6 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_nonsymmetric>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
# TODO use symmetric memory above if doing regional MOM6
# Note that the explicit dependency on OpenMP is only because FMS is not declaring its OpenMP dependency correctly
# This will be fixed with fms version 2024.01
target_link_libraries(mom6 PUBLIC NetCDF::NetCDF_C NetCDF::NetCDF_Fortran fms OpenMP::OpenMP_Fortran OpenMP::OpenMP_C)

####################################################
### Create standalone MOM6 executable
####################################################
if(MOM6SOLO)
  message("Building MOM6 standalone executable")
  add_executable(mom6solo ${mom6_solo_src_files})
  add_dependencies(mom6solo mom6)
  set_target_properties(mom6solo PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_solo)
  set_target_properties(mom6solo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  target_include_directories(mom6solo PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
  target_include_directories(mom6solo PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_solo>)
  target_include_directories(mom6solo PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_nonsymmetric>
                                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
  target_link_libraries(mom6solo PUBLIC
    NetCDF::NetCDF_C
    NetCDF::NetCDF_Fortran
    fms
    mom6 )

endif()

####################################################
### Install
####################################################
install (
    TARGETS mom6
    EXPORT mom6-config
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    COMPONENT Library)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/MOM6)
install(EXPORT      mom6-config
        DESTINATION lib/cmake)

