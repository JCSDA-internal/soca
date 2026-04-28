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

std::string layerName(const std::string & base, std::size_t l) {
  // CICE convention: qice001..qice00N, sice001..sice00N, qsno001..qsno00N.
  std::string out = base;
  if (l < 10) {
    out += "00";
    out += static_cast<char>('0' + static_cast<int>(l));
  } else if (l < 100) {
    out += "0";
    out += static_cast<char>('0' + static_cast<int>(l / 10));
    out += static_cast<char>('0' + static_cast<int>(l % 10));
  } else {
    // Should not happen for ice or snow.
    out += std::to_string(l);
  }
  return out;
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

// ---------------------------------------------------------------------------
// ownedGlobalIndices
// ---------------------------------------------------------------------------
std::vector<atlas::gidx_t> CiceRestartIO::ownedGlobalIndices() const {
  const auto & fs    = geom_.functionSpace();
  const auto ghost   = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindex  = atlas::array::make_view<atlas::gidx_t, 1>(
                          fs.global_index());
  std::vector<atlas::gidx_t> owned;
  owned.reserve(ghost.size());
  for (std::size_t jnode = 0; jnode < ghost.size(); ++jnode) {
    if (ghost(jnode)) continue;
    const atlas::gidx_t gidx = gindex(jnode);
    if (gidx > 0) owned.push_back(gidx);
  }
  return owned;
}

// ---------------------------------------------------------------------------
// readThermo
//
// Strategy: each rank announces its owned-gidx list to root (gather-of-int
// counts + gather-of-gidx data). Root reads each thermo field globally from
// the input file, then for each rank slices out the per-node values in the
// rank's order and sends them back via a single scatterv. Each variable is
// scattered independently.
// ---------------------------------------------------------------------------
CiceRestartIO::ThermoFrame CiceRestartIO::readThermo(
    std::size_t ncat, std::size_t iceLev, std::size_t snoLev) const {
  const auto & comm = geom_.getComm();
  const std::size_t nr = comm.size();

  ThermoFrame frame;
  frame.ncat   = ncat;
  frame.iceLev = iceLev;
  frame.snoLev = snoLev;

  // Local owned gidx list and count.
  std::vector<atlas::gidx_t> ownedGidx = ownedGlobalIndices();
  const int localCount = static_cast<int>(ownedGidx.size());
  frame.nOwnedNodes = ownedGidx.size();

  // Allocate the local buffers.
  frame.Tsfcn.assign(frame.nOwnedNodes * ncat,         0.0);
  frame.qice .assign(frame.nOwnedNodes * ncat * iceLev, 0.0);
  frame.sice .assign(frame.nOwnedNodes * ncat * iceLev, 0.0);
  frame.qsno .assign(frame.nOwnedNodes * ncat * snoLev, 0.0);
  frame.apnd .assign(frame.nOwnedNodes * ncat,         0.0);
  frame.hpnd .assign(frame.nOwnedNodes * ncat,         0.0);
  frame.ipnd .assign(frame.nOwnedNodes * ncat,         0.0);

  // Gather per-rank counts then per-rank gidx lists onto root.
  std::vector<int> counts(nr, 0);
  comm.allGather(localCount, counts.begin(), counts.end());
  std::vector<int> displs(nr, 0);
  std::size_t totalNodes = 0;
  for (std::size_t r = 0; r < nr; ++r) {
    displs[r] = static_cast<int>(totalNodes);
    totalNodes += static_cast<std::size_t>(counts[r]);
  }
  // Convert local gidx to int for transport (CICE grid fits in int easily).
  std::vector<int> sendGidx(localCount);
  for (int i = 0; i < localCount; ++i) {
    sendGidx[i] = static_cast<int>(ownedGidx[i]);
  }
  std::vector<int> allGidx;
  if (comm.rank() == 0) allGidx.resize(totalNodes);
  comm.gatherv(sendGidx, allGidx, counts, displs, 0);

  // Open the input restart on root.
  int ncid = -1;
  std::size_t ni = 0, nj = 0;
  if (comm.rank() == 0) {
    std::size_t ncatChk = 0;
    readDims(ncatChk, nj, ni);
    if (ncatChk != ncat) {
      throw eckit::UserError(
          "CiceRestartIO::readThermo: ncat mismatch", Here());
    }
    ncCheck(nc_open(inputFile_.c_str(), NC_NOWRITE, &ncid),
            "nc_open input thermo");
  } else {
    // Other ranks need ni, nj for the broadcast pattern below.
    std::size_t ncatChk = 0;
    readDims(ncatChk, nj, ni);
  }

  auto scatterCatField = [&](const std::string & varName,
                             std::vector<double> & local) {
    // Read globally on root.
    std::vector<double> global;
    if (comm.rank() == 0) {
      global.assign(ncat * nj * ni, 0.0);
      const int varid = inqVarOptional(ncid, varName);
      if (varid >= 0) {
        ncCheck(nc_get_var_double(ncid, varid, global.data()),
                ("nc_get_var " + varName).c_str());
      }
    }
    // Build send buffer per rank: for each rank r, gather its nodes' values
    // (k-major, so each node contributes ncat doubles).
    const std::size_t recSize = ncat;
    std::vector<double> send;
    if (comm.rank() == 0) {
      send.resize(totalNodes * recSize, 0.0);
      for (std::size_t r = 0; r < nr; ++r) {
        for (int idx = 0; idx < counts[r]; ++idx) {
          const int gidx = allGidx[displs[r] + idx];
          if (gidx <= 0) continue;
          const std::size_t lin = static_cast<std::size_t>(gidx - 1);
          const std::size_t i = lin % ni;
          const std::size_t j = lin / ni;
          if (j >= nj) continue;
          for (std::size_t k = 0; k < ncat; ++k) {
            send[(static_cast<std::size_t>(displs[r]) + idx) * recSize + k] =
                global[(k * nj + j) * ni + i];
          }
        }
      }
    }
    // Per-rank send/recv counts in doubles.
    std::vector<int> sCounts(nr), sDispls(nr);
    for (std::size_t r = 0; r < nr; ++r) {
      sCounts[r] = counts[r] * static_cast<int>(recSize);
      sDispls[r] = displs[r] * static_cast<int>(recSize);
    }
    comm.scatterv(send.begin(), send.end(),
                  sCounts, sDispls,
                  local.begin(), local.end(), 0);
  };

  auto scatterCatLevField = [&](const std::string & baseName,
                                std::size_t nlyr,
                                std::vector<double> & local) {
    const std::size_t recSize = ncat * nlyr;
    std::vector<double> send;
    if (comm.rank() == 0) {
      send.assign(totalNodes * recSize, 0.0);
    }
    // Read each layer in turn on root, fold its values into the send buffer.
    for (std::size_t l = 0; l < nlyr; ++l) {
      std::vector<double> global;
      if (comm.rank() == 0) {
        global.assign(ncat * nj * ni, 0.0);
        const int varid = inqVarOptional(ncid, layerName(baseName, l + 1));
        if (varid >= 0) {
          ncCheck(nc_get_var_double(ncid, varid, global.data()),
                  ("nc_get_var " + layerName(baseName, l + 1)).c_str());
        }
        for (std::size_t r = 0; r < nr; ++r) {
          for (int idx = 0; idx < counts[r]; ++idx) {
            const int gidx = allGidx[displs[r] + idx];
            if (gidx <= 0) continue;
            const std::size_t lin = static_cast<std::size_t>(gidx - 1);
            const std::size_t i = lin % ni;
            const std::size_t j = lin / ni;
            if (j >= nj) continue;
            for (std::size_t k = 0; k < ncat; ++k) {
              send[(static_cast<std::size_t>(displs[r]) + idx) * recSize
                   + k * nlyr + l] = global[(k * nj + j) * ni + i];
            }
          }
        }
      }
    }
    std::vector<int> sCounts(nr), sDispls(nr);
    for (std::size_t r = 0; r < nr; ++r) {
      sCounts[r] = counts[r] * static_cast<int>(recSize);
      sDispls[r] = displs[r] * static_cast<int>(recSize);
    }
    comm.scatterv(send.begin(), send.end(),
                  sCounts, sDispls,
                  local.begin(), local.end(), 0);
  };

  scatterCatField("Tsfcn", frame.Tsfcn);
  scatterCatLevField("qice", iceLev, frame.qice);
  scatterCatLevField("sice", iceLev, frame.sice);
  scatterCatLevField("qsno", snoLev, frame.qsno);
  scatterCatField("apnd", frame.apnd);
  scatterCatField("hpnd", frame.hpnd);
  scatterCatField("ipnd", frame.ipnd);

  if (comm.rank() == 0) {
    ncCheck(nc_close(ncid), "nc_close input thermo");
  }
  comm.barrier();
  return frame;
}

// ---------------------------------------------------------------------------
// flushThermo
//
// Mirror of readThermo: each rank packs its (gidx, [vals]) records, gatherv
// to root, root reassembles globally and overwrites the variables in the
// already-existing output file.
// ---------------------------------------------------------------------------
void CiceRestartIO::flushThermo(const ThermoFrame & frame) const {
  const auto & comm = geom_.getComm();
  const std::size_t nr = comm.size();

  std::vector<atlas::gidx_t> ownedGidx = ownedGlobalIndices();
  const int localCount = static_cast<int>(ownedGidx.size());
  if (frame.nOwnedNodes != ownedGidx.size()) {
    throw eckit::UserError(
        "CiceRestartIO::flushThermo: frame mismatch", Here());
  }

  std::vector<int> counts(nr, 0);
  comm.allGather(localCount, counts.begin(), counts.end());
  std::vector<int> displs(nr, 0);
  std::size_t totalNodes = 0;
  for (std::size_t r = 0; r < nr; ++r) {
    displs[r] = static_cast<int>(totalNodes);
    totalNodes += static_cast<std::size_t>(counts[r]);
  }
  std::vector<int> sendGidx(localCount);
  for (int i = 0; i < localCount; ++i) {
    sendGidx[i] = static_cast<int>(ownedGidx[i]);
  }
  std::vector<int> allGidx;
  if (comm.rank() == 0) allGidx.resize(totalNodes);
  comm.gatherv(sendGidx, allGidx, counts, displs, 0);

  std::size_t ncatChk = 0, ni = 0, nj = 0;
  readDims(ncatChk, nj, ni);

  // Open the output file on root for write-overwrite.
  int ncid = -1;
  if (comm.rank() == 0) {
    ncCheck(nc_open(outputFile_.c_str(), NC_WRITE, &ncid),
            "nc_open output thermo");
  }

  auto gatherAndWriteCat = [&](const std::string & varName,
                               const std::vector<double> & local) {
    const std::size_t recSize = frame.ncat;
    std::vector<int> sCounts(nr), sDispls(nr);
    for (std::size_t r = 0; r < nr; ++r) {
      sCounts[r] = counts[r] * static_cast<int>(recSize);
      sDispls[r] = displs[r] * static_cast<int>(recSize);
    }
    std::vector<double> recv;
    if (comm.rank() == 0) recv.resize(totalNodes * recSize, 0.0);
    comm.gatherv(local, recv, sCounts, sDispls, 0);
    if (comm.rank() == 0) {
      const int varid = inqVarOptional(ncid, varName);
      if (varid < 0) return;
      std::vector<double> global(frame.ncat * nj * ni, 0.0);
      // Default-fill global from the existing file so any unowned cell
      // (land outside the partition, rounding) keeps its template value.
      ncCheck(nc_get_var_double(ncid, varid, global.data()),
              ("nc_get_var (template) " + varName).c_str());
      for (std::size_t r = 0; r < nr; ++r) {
        for (int idx = 0; idx < counts[r]; ++idx) {
          const int gidx = allGidx[displs[r] + idx];
          if (gidx <= 0) continue;
          const std::size_t lin = static_cast<std::size_t>(gidx - 1);
          const std::size_t i = lin % ni;
          const std::size_t j = lin / ni;
          if (j >= nj) continue;
          for (std::size_t k = 0; k < frame.ncat; ++k) {
            global[(k * nj + j) * ni + i] = recv[
                (static_cast<std::size_t>(displs[r]) + idx) * recSize + k];
          }
        }
      }
      ncCheck(nc_put_var_double(ncid, varid, global.data()),
              ("nc_put_var " + varName).c_str());
    }
  };

  auto gatherAndWriteCatLev = [&](const std::string & baseName,
                                  std::size_t nlyr,
                                  const std::vector<double> & local) {
    const std::size_t recSize = frame.ncat * nlyr;
    std::vector<int> sCounts(nr), sDispls(nr);
    for (std::size_t r = 0; r < nr; ++r) {
      sCounts[r] = counts[r] * static_cast<int>(recSize);
      sDispls[r] = displs[r] * static_cast<int>(recSize);
    }
    std::vector<double> recv;
    if (comm.rank() == 0) recv.resize(totalNodes * recSize, 0.0);
    comm.gatherv(local, recv, sCounts, sDispls, 0);
    if (comm.rank() == 0) {
      for (std::size_t l = 0; l < nlyr; ++l) {
        const std::string varName = layerName(baseName, l + 1);
        const int varid = inqVarOptional(ncid, varName);
        if (varid < 0) continue;
        std::vector<double> global(frame.ncat * nj * ni, 0.0);
        ncCheck(nc_get_var_double(ncid, varid, global.data()),
                ("nc_get_var (template) " + varName).c_str());
        for (std::size_t r = 0; r < nr; ++r) {
          for (int idx = 0; idx < counts[r]; ++idx) {
            const int gidx = allGidx[displs[r] + idx];
            if (gidx <= 0) continue;
            const std::size_t lin = static_cast<std::size_t>(gidx - 1);
            const std::size_t i = lin % ni;
            const std::size_t j = lin / ni;
            if (j >= nj) continue;
            for (std::size_t k = 0; k < frame.ncat; ++k) {
              global[(k * nj + j) * ni + i] = recv[
                  (static_cast<std::size_t>(displs[r]) + idx) * recSize
                  + k * nlyr + l];
            }
          }
        }
        ncCheck(nc_put_var_double(ncid, varid, global.data()),
                ("nc_put_var " + varName).c_str());
      }
    }
  };

  gatherAndWriteCat("Tsfcn", frame.Tsfcn);
  gatherAndWriteCatLev("qice", frame.iceLev, frame.qice);
  gatherAndWriteCatLev("sice", frame.iceLev, frame.sice);
  gatherAndWriteCatLev("qsno", frame.snoLev, frame.qsno);
  gatherAndWriteCat("apnd", frame.apnd);
  gatherAndWriteCat("hpnd", frame.hpnd);
  gatherAndWriteCat("ipnd", frame.ipnd);

  if (comm.rank() == 0) {
    ncCheck(nc_close(ncid), "nc_close output thermo");
  }
  comm.barrier();
}

}  // namespace soca
