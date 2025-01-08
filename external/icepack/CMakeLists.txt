# The file was modified from the original wrapper in the ufs-weather-model
# https://github.com/ufs-community/ufs-weather-model/blob/develop/CICE-interface/CMakeLists.txt

# Too many files to list, so include them via this file
include("icepack_files.cmake")

# Collect source files for library
list(APPEND lib_src_files
  ${icepack_files})

### Create target library and set PUBLIC interfaces on the library
add_library(icepack STATIC ${lib_src_files})
set_target_properties(icepack PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
#target_compile_definitions(icepack PUBLIC "${_cice_defs}")
target_include_directories(icepack PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)

###############################################################################
### Install
###############################################################################

install(
  TARGETS icepack
  EXPORT  icepack-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/icepack)

install(EXPORT      icepack-config
        DESTINATION lib/cmake)