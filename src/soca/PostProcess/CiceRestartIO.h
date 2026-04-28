/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "eckit/mpi/Comm.h"

namespace soca {

class Geometry;

// ---------------------------------------------------------------------------
/// @brief Write a CICE-format sea-ice restart NetCDF from per-category atlas
/// fields on the soca/MOM6 grid.
///
/// The CICE restart layout is `aicen[ncat,nj,ni]` (and similar for vicen,
/// vsnon, Tsfcn) plus per-layer `qice00N[ncat,nj,ni]`, `sice00N[ncat,nj,ni]`,
/// `qsno00N[ncat,nj,ni]`. To preserve all the unmodified CICE fields and
/// global attributes (apnd, hpnd, ipnd, iceumask, istep1, ...) we copy the
/// input restart file to the output path and overwrite only the variables
/// for which we have updated fields in the FieldSet.
///
/// Parallel layout: the soca grid is the same (ni x nj) as CICE, distributed
/// across MPI ranks via MOM6's domain decomposition. Each rank holds owned
/// nodes whose global indices map back to (i,j) on the global grid. We
/// gather (i, j, value) triples to root and assemble the global field there.
// ---------------------------------------------------------------------------
class CiceRestartIO {
 public:
  CiceRestartIO(const Geometry & geom,
                const std::string & inputFile,
                const std::string & outputFile);

  /// Write a CICE restart by overwriting variables in a copy of the input
  /// restart. `fset` must contain the per-category fields named
  /// `sea_ice_categoryK_area_fraction`, `sea_ice_categoryK_volume`,
  /// `sea_ice_snow_categoryK_volume` for K = 1..ncat.
  void write(const atlas::FieldSet & fset, std::size_t ncat) const;

 private:
  const Geometry & geom_;
  std::string inputFile_;
  std::string outputFile_;

  /// Gather a single per-category field from local FieldSet nodes onto rank 0
  /// in CICE-layout order (k-major, then j*ni + i). Returns the global array
  /// on rank 0; empty vector on other ranks.
  std::vector<double> gatherCategoryField(
      const std::vector<atlas::array::ArrayView<double, 2>> & catViews,
      std::size_t ni, std::size_t nj, std::size_t ncat) const;

  /// Discover (ncat, nj, ni) from the input file (root reads, then bcasts).
  void readDims(std::size_t & ncat, std::size_t & nj, std::size_t & ni) const;
};

}  // namespace soca
