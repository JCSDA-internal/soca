/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/PostProcess/CiceRestartIO.h"

#include <netcdf.h>

#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"

namespace soca {

namespace {

// Throw on NetCDF error.
void ncCheck(int rc, const char * where) {
  if (rc != NC_NOERR) {
    std::ostringstream oss;
    oss << "CiceRestartIO: " << where << ": " << nc_strerror(rc);
    throw eckit::Exception(oss.str(), Here());
  }
}

// Cheap binary file copy. NetCDF formats are self-contained, so a byte copy
// of the input file gives a writable template for nf90_open(NC_WRITE).
void copyFile(const std::string & src, const std::string & dst) {
  std::ifstream in(src, std::ios::binary);
  if (!in) {
    throw eckit::CantOpenFile(src);
  }
  std::ofstream out(dst, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw eckit::CantOpenFile(dst);
  }
  out << in.rdbuf();
}

// Look up a variable; return its varid, or -1 if it doesn't exist.
int inqVarOptional(int ncid, const std::string & name) {
  int varid = -1;
  const int rc = nc_inq_varid(ncid, name.c_str(), &varid);
  if (rc == NC_ENOTVAR) return -1;
  ncCheck(rc, ("nc_inq_varid " + name).c_str());
  return varid;
}

}  // namespace

// ---------------------------------------------------------------------------

CiceRestartIO::CiceRestartIO(const Geometry & geom,
                             const std::string & inputFile,
                             const std::string & outputFile)
  : geom_(geom), inputFile_(inputFile), outputFile_(outputFile) {}

// ---------------------------------------------------------------------------

void CiceRestartIO::readDims(std::size_t & ncat,
                             std::size_t & nj,
                             std::size_t & ni) const {
  const auto & comm = geom_.getComm();
  std::array<std::size_t, 3> dims{0, 0, 0};
  if (comm.rank() == 0) {
    int ncid = -1;
    ncCheck(nc_open(inputFile_.c_str(), NC_NOWRITE, &ncid), "nc_open input");
    int dimid = -1;
    std::size_t len = 0;
    ncCheck(nc_inq_dimid(ncid, "ncat", &dimid), "inq ncat");
    ncCheck(nc_inq_dimlen(ncid, dimid, &len), "len ncat");
    dims[0] = len;
    ncCheck(nc_inq_dimid(ncid, "nj", &dimid), "inq nj");
    ncCheck(nc_inq_dimlen(ncid, dimid, &len), "len nj");
    dims[1] = len;
    ncCheck(nc_inq_dimid(ncid, "ni", &dimid), "inq ni");
    ncCheck(nc_inq_dimlen(ncid, dimid, &len), "len ni");
    dims[2] = len;
    ncCheck(nc_close(ncid), "nc_close input");
  }
  comm.broadcast(dims.begin(), dims.end(), 0);
  ncat = dims[0];
  nj   = dims[1];
  ni   = dims[2];
}

// ---------------------------------------------------------------------------

std::vector<double> CiceRestartIO::gatherCategoryField(
    const std::vector<atlas::array::ArrayView<double, 2>> & catViews,
    std::size_t ni, std::size_t nj, std::size_t ncat) const {
  const auto & comm = geom_.getComm();
  const auto & fs   = geom_.functionSpace();

  // Local owned (i, j) and per-cat values.
  const auto ghost  = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindex = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());
  const std::size_t nnodes = ghost.size();

  // Pack: each owned node -> [gidx-1] index, then ncat values.
  // Use a single contiguous double buffer with the gidx encoded as the first
  // entry of each record. (gidx fits in double for any practical CICE grid.)
  const std::size_t recSize = 1 + ncat;
  std::vector<double> sendbuf;
  sendbuf.reserve(recSize * nnodes);
  for (std::size_t jnode = 0; jnode < nnodes; ++jnode) {
    if (ghost(jnode)) continue;
    const atlas::gidx_t gidx = gindex(jnode);
    if (gidx <= 0) continue;
    sendbuf.push_back(static_cast<double>(gidx));
    for (std::size_t k = 0; k < ncat; ++k) {
      sendbuf.push_back(catViews[k](jnode, 0));
    }
  }
  const int sendCount = static_cast<int>(sendbuf.size());

  // Gather counts then values to root.
  std::vector<int> counts(comm.size(), 0);
  comm.allGather(sendCount, counts.begin(), counts.end());
  std::vector<int> displs(comm.size(), 0);
  std::size_t total = 0;
  for (std::size_t r = 0; r < comm.size(); ++r) {
    displs[r] = static_cast<int>(total);
    total += static_cast<std::size_t>(counts[r]);
  }
  std::vector<double> recvbuf;
  if (comm.rank() == 0) recvbuf.resize(total);
  comm.gatherv(sendbuf, recvbuf, counts, displs, 0);

  if (comm.rank() != 0) return {};

  // Assemble (ncat, nj, ni) global array in CICE layout (k-major, then j, i).
  std::vector<double> global(ncat * nj * ni, 0.0);
  for (std::size_t off = 0; off + recSize <= recvbuf.size(); off += recSize) {
    const auto gidx = static_cast<std::int64_t>(recvbuf[off]);
    if (gidx <= 0) continue;
    const std::size_t lin = static_cast<std::size_t>(gidx - 1);  // 0-based
    const std::size_t i   = lin % ni;
    const std::size_t j   = lin / ni;
    if (j >= nj) continue;
    for (std::size_t k = 0; k < ncat; ++k) {
      global[(k * nj + j) * ni + i] = recvbuf[off + 1 + k];
    }
  }
  return global;
}

// ---------------------------------------------------------------------------

void CiceRestartIO::write(const atlas::FieldSet & fset, std::size_t ncatIn) const {
  std::size_t ncat = 0, nj = 0, ni = 0;
  readDims(ncat, nj, ni);
  if (ncat != ncatIn) {
    std::ostringstream oss;
    oss << "CiceRestartIO::write: ncat in restart (" << ncat
        << ") != ncat from postprocess (" << ncatIn << ")";
    throw eckit::UserError(oss.str(), Here());
  }

  const auto & comm = geom_.getComm();

  // Copy the input restart to the output path on root, then broadcast that
  // it's ready. Other ranks don't touch the output file.
  if (comm.rank() == 0) {
    copyFile(inputFile_, outputFile_);
  }
  comm.barrier();

  // For each tracked variable group, gather its per-cat field and write.
  // Only `aicen`, `vicen`, `vsnon` are produced by Stage A/B today; thermo
  // (Tsfcn, qiceNN, siceNN, qsnoNN) and pond (apnd/hpnd/ipnd) fields will be
  // added once Stage C is wired.
  const std::vector<std::pair<std::string, std::string>> groups{
    {"aicen", "sea_ice_category{}_area_fraction"},
    {"vicen", "sea_ice_category{}_volume"},
    {"vsnon", "sea_ice_snow_category{}_volume"},
  };

  for (const auto & [ncVarName, fieldPattern] : groups) {
    std::vector<atlas::array::ArrayView<double, 2>> catViews;
    catViews.reserve(ncat);
    for (std::size_t k = 0; k < ncat; ++k) {
      const std::string fieldName = std::string(fieldPattern).replace(
          fieldPattern.find("{}"), 2, std::to_string(k + 1));
      catViews.push_back(
          atlas::array::make_view<double, 2>(fset.field(fieldName)));
    }
    auto global = gatherCategoryField(catViews, ni, nj, ncat);

    if (comm.rank() == 0) {
      int ncid = -1;
      ncCheck(nc_open(outputFile_.c_str(), NC_WRITE, &ncid),
              "nc_open output");
      const int varid = inqVarOptional(ncid, ncVarName);
      if (varid < 0) {
        oops::Log::warning() << "CiceRestartIO: variable '" << ncVarName
                             << "' not found in restart; skipping" << std::endl;
      } else {
        ncCheck(nc_put_var_double(ncid, varid, global.data()),
                ("nc_put_var " + ncVarName).c_str());
      }
      ncCheck(nc_close(ncid), "nc_close output");
    }
  }

  comm.barrier();
}

}  // namespace soca
