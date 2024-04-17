/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include "atlas/interpolation/Interpolation.h"
#include "atlas/util/Point.h"

#include "soca/Utils/readNcAndInterp.h"


namespace soca {

atlas::FieldSet readNcAndInterp(
    const std::string & filename,
    const std::vector<std::string> & vars,
    const atlas::FunctionSpace & dstFunctionSpace ) {
  atlas::FieldSet fieldSet;

  int ncid;
  if (nc_open(filename.c_str(), NC_NOWRITE, &ncid)) {
    throw std::runtime_error("Failed to open NetCDF file.");
  }

  // Read latitude and longitude
  int latVarId, lonVarId, nDims, dimId;
  nc_inq_varid(ncid, "latitude", &latVarId);
  nc_inq_varid(ncid, "longitude", &lonVarId);

  size_t srcSize;
  nc_inq_varndims(ncid, latVarId, &nDims);
  ASSERT(nDims == 1);
  nc_inq_vardimid(ncid, latVarId, &dimId);
  nc_inq_dimlen(ncid, dimId, &srcSize);

  std::vector<double> latitudes(srcSize), longitudes(srcSize);
  nc_get_var_double(ncid, latVarId, latitudes.data());
  nc_get_var_double(ncid, lonVarId, longitudes.data());

  // create src functionspace
  auto srcLonLatField = atlas::Field("lonlat",
    atlas::array::make_datatype<double>(), atlas::array::make_shape(srcSize, 2));
  auto srcLonLatView = atlas::array::make_view<double, 2>(srcLonLatField);
  for (size_t i = 0; i < srcSize; i++) {
    auto point = atlas::PointLonLat(longitudes[i], latitudes[i]);
    point.normalise();
    srcLonLatView(i, 0) = point.lon();
    srcLonLatView(i, 1) = point.lat();
  }
  const auto srcFunctionSpace = atlas::functionspace::PointCloud(srcLonLatField);

  // create interpolation
  eckit::LocalConfiguration interpConfig;
  interpConfig.set("type", "k-nearest-neighbours");
  interpConfig.set("k-nearest-neighbours", 10); // do we really need this many?
  atlas::Interpolation interp(interpConfig, srcFunctionSpace, dstFunctionSpace);

  // read and interpolate fields
  for (const std::string & varName : vars) {
    // Read variable data
    int varId;
    nc_inq_varid(ncid, varName.c_str(), &varId);
    std::vector<double> varData(srcSize);
    nc_get_var_double(ncid, varId, varData.data());

    // Create a field for the variable
    atlas::Field sourceField = srcFunctionSpace.createField<double>(
    atlas::option::name(varName) | atlas::option::levels(1));
    auto view = atlas::array::make_view<double, 2>(sourceField);
    for (size_t i = 0; i < varData.size(); ++i) {
      view(i, 0) = varData[i];
    }

    // interpolate
    atlas::Field dstField = dstFunctionSpace.createField<double>(
      atlas::option::name(varName) | atlas::option::levels(1));
    interp.execute(sourceField, dstField);
    fieldSet.add(dstField);
  }
  return fieldSet;
}

}  // namespace soca
