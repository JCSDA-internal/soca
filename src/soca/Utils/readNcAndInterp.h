/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#pragma once

#include <string>
#include <vector>

#include "atlas/functionspace.h"

namespace soca
{

/// @brief Read a netCDF file and interpolate the values to the destination atlas FunctionSpace.
///
/// It is assumed that variables in the input file are all 1D, and that "latitude" and "longitude"
/// variables are present in addition to the variables listed in `vars`. All  resulting fields are 2D
///
/// @param filename The name of the netCDF file to read.
/// @param vars A list of variable names to read in and process
/// @param dstFunctionSpace The destination atlas FunctionSpace that input variables are
///        interpolated to.
/// @return An atlas FieldSet with the interpolated fields.
atlas::FieldSet readNcAndInterp(
  const std::string & filename,
  const std::vector<std::string> & vars,
  const atlas::FunctionSpace & dstFunctionSpace);

}  // namespace soca

