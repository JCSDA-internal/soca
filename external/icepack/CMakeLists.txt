# The file was modified from the original wrapper in the ufs-weather-model
# https://github.com/ufs-community/ufs-weather-model/blob/develop/CICE-interface/CMakeLists.txt

# Too many files to list, so include them via this file
include("icepack_files.cmake")

# Collect source files for library
list(APPEND lib_src_files
  ${icepack_files})

# Define a variable for the library name
# (NOTE: we do not call this icepack, because, for now, only a subset of it
#  is being built based on what soca needs)
set(LIBRARY_NAME icepack_limited)

## Create target library and set PUBLIC interfaces on the library
add_library(${LIBRARY_NAME} STATIC ${lib_src_files})
set_target_properties(${LIBRARY_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

target_include_directories(${LIBRARY_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
  $<INSTALL_INTERFACE:mod>)

# ##############################################################################
# ## Install
# ##############################################################################
install(
  TARGETS ${LIBRARY_NAME}
  EXPORT ${LIBRARY_NAME}-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBRARY_NAME})

install(EXPORT ${LIBRARY_NAME}-config
  DESTINATION lib/cmake)