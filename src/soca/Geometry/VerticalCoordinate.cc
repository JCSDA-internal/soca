/*
 * (C) Copyright 2026- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/missingValues.h"

#include "soca/Geometry/VerticalCoordinate.h"

namespace soca {

namespace {

bool wet(const double h) {return h != util::missingValue<double>() && h > 0.0;}

double median(std::vector<double> & v) {
  const size_t m = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + m, v.end());
  if (v.size() % 2 == 1) return v[m];
  return 0.5 * (*std::max_element(v.begin(), v.begin() + m) + v[m]);
}

void ncCheck(const int status, const std::string & what) {
  if (status != NC_NOERR) {
    throw eckit::UserError(what + ": " + nc_strerror(status), Here());
  }
}

}  // namespace

// -----------------------------------------------------------------------------
std::vector<double> verticalCoordinate(const std::string & units, const size_t nlevs,
                                       const std::vector<double> & layerDepths,
                                       const double logDepthOffset) {
  if (units == "levels") {
    std::vector<double> coords(nlevs);
    std::iota(coords.begin(), coords.end(), 1.0);
    return coords;
  }

  if (units != "m" && units != "log depth") {
    throw eckit::NotImplemented("Vertical coordinate '" + units + "' not supported. "
      "Supported: 'levels', 'm', 'log depth'.", Here());
  }

  if (layerDepths.size() != nlevs) {
    throw eckit::UserError("Vertical coordinate '" + units + "' needs the layer depths "
      "('layer thickness file' in the geometry) with one value per layer ("
      + std::to_string(nlevs) + "), got " + std::to_string(layerDepths.size()), Here());
  }
  for (size_t k = 0; k < nlevs; ++k) {
    if (layerDepths[k] < 0.0 || (k > 0 && layerDepths[k] <= layerDepths[k - 1])) {
      throw eckit::UserError("The layer depths must be non-negative and increasing", Here());
    }
  }
  if (units == "m") return layerDepths;

  if (logDepthOffset <= 0.0) {
    throw eckit::UserError("'log depth offset' must be positive", Here());
  }
  std::vector<double> coords(nlevs);
  for (size_t k = 0; k < nlevs; ++k) coords[k] = std::log(layerDepths[k] + logDepthOffset);
  return coords;
}

// -----------------------------------------------------------------------------
namespace {

std::vector<double> readLayerDepths(const std::string & file, const std::string & var,
                                    const size_t nlevs) {
  int ncid;
  ncCheck(nc_open(file.c_str(), NC_NOWRITE, &ncid), "Cannot open '" + file + "'");
  struct Closer {
    int id;
    ~Closer() {nc_close(id);}
  } closer{ncid};

  int varid, ndims;
  ncCheck(nc_inq_varid(ncid, var.c_str(), &varid), "No variable '" + var + "' in '" + file + "'");
  ncCheck(nc_inq_varndims(ncid, varid, &ndims), "'" + var + "'");
  if (ndims != 3 && ndims != 4) {
    throw eckit::UserError("'" + var + "' in '" + file + "' must have dimensions "
                           "([time,] layer, y, x)", Here());
  }
  std::vector<int> dimids(ndims);
  ncCheck(nc_inq_vardimid(ncid, varid, dimids.data()), "'" + var + "'");
  std::vector<size_t> count(ndims);
  for (int i = 0; i < ndims; ++i) {
    ncCheck(nc_inq_dimlen(ncid, dimids[i], &count[i]), "'" + var + "'");
  }
  const int kz = ndims - 3;
  if (count[kz] != nlevs) {
    throw eckit::UserError("'" + var + "' in '" + file + "' has " + std::to_string(count[kz])
      + " layers, the geometry " + std::to_string(nlevs), Here());
  }
  if (ndims == 4 && count[0] == 0) {
    throw eckit::UserError("'" + var + "' in '" + file + "' has no time record", Here());
  }

  // values marked missing in the file: _FillValue (or the netCDF default) and missing_value
  std::vector<double> fills;
  double att;
  if (nc_get_att_double(ncid, varid, "_FillValue", &att) == NC_NOERR) {
    fills.push_back(att);
  } else {
    nc_type type;
    ncCheck(nc_inq_vartype(ncid, varid, &type), "'" + var + "'");
    fills.push_back(type == NC_FLOAT ? static_cast<double>(NC_FILL_FLOAT) : NC_FILL_DOUBLE);
  }
  if (nc_get_att_double(ncid, varid, "missing_value", &att) == NC_NOERR) fills.push_back(att);

  std::vector<size_t> start(ndims, 0);
  if (ndims == 4) count[0] = 1;
  count[kz] = 1;
  const size_t ncol = count[kz + 1] * count[kz + 2];
  auto readLayer = [&](const size_t k, std::vector<double> & h) {
    start[kz] = k;
    ncCheck(nc_get_vara_double(ncid, varid, start.data(), count.data(), h.data()),
            "Reading '" + var + "' from '" + file + "'");
    for (double & v : h) {
      if (std::find(fills.begin(), fills.end(), v) != fills.end()) {
        v = util::missingValue<double>();
      }
    }
  };

  // layer k needs layer k+1 to leave out partial cells at the sea floor
  std::vector<double> depths(nlevs), top(ncol, 0.0), h(ncol), below(ncol), centres;
  readLayer(0, h);
  for (size_t k = 0; k < nlevs; ++k) {
    const bool bottom = (k + 1 == nlevs);
    if (!bottom) readLayer(k + 1, below);
    centres.clear();
    for (size_t i = 0; i < ncol; ++i) {
      if (wet(h[i]) && (bottom || wet(below[i]))) centres.push_back(top[i] + 0.5 * h[i]);
    }
    if (centres.empty()) {
      throw eckit::UserError("'" + var + "' in '" + file + "': layer " + std::to_string(k + 1)
                             + " has no wet column", Here());
    }
    depths[k] = median(centres);
    for (size_t i = 0; i < ncol; ++i) {
      if (wet(h[i])) top[i] += h[i];
    }
    std::swap(h, below);
  }
  return depths;
}

}  // namespace

// -----------------------------------------------------------------------------
std::vector<double> layerDepthsFromFile(const std::string & file, const std::string & var,
                                        const size_t nlevs, const eckit::mpi::Comm & comm) {
  std::vector<double> depths(nlevs, 0.0);
  std::string error;
  if (comm.rank() == 0) {
    try {
      depths = readLayerDepths(file, var, nlevs);
    } catch (const std::exception & e) {
      error = e.what();
    }
  }
  int failed = error.empty() ? 0 : 1;
  comm.broadcast(failed, 0);
  if (failed) {
    throw eckit::UserError("Layer depths from '" + file + "' failed"
                           + (error.empty() ? std::string(" on rank 0") : ": " + error), Here());
  }
  comm.broadcast(depths.data(), nlevs, 0);
  return depths;
}

// -----------------------------------------------------------------------------

}  // namespace soca
