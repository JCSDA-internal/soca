# (C) Copyright 2026 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


# Set compiler flags for basic build types,
# for compilers where this is not provided by ecbuild.
include(build_type_compiler_flags)

# Set JEDI's common compiler flags
include(jedi_common_compiler_flags)

# Set SOCA-specific compiler flags
add_definitions(-Duse_libMPI -Duse_netCDF -DSPMD)
if(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
  add_definitions( -DC_F_NO_DEALLOC )
endif()
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  ecbuild_add_fortran_flags("-ffree-line-length-none")
  ecbuild_add_fortran_flags("-ffpe-trap=invalid,zero,overflow,underflow" BUILD DEBUG)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)  # Intel or IntelLLVM
  # 1. icepack raises FPEs when compiled with default FP settings => fixable with =precise
  # 2. some soca ctest output needs =strict to remain consistent between Intel and IntelLLVM
  #    TODO: go back to =precise and adjust test tolerances?
  ecbuild_add_fortran_flags("-fp-model=strict")
endif()
