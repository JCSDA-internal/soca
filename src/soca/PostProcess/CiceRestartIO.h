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

  // -------------------------------------------------------------------------
  /// @brief Per-rank container for the CICE restart variables that don't
  /// flow through the soca State / atlas FieldSet (thermo + ponds).
  ///
  /// Each variable is held as a vector of doubles indexed by
  /// `node * ncat + k` (or `node * ncat * nlyr + k * nlyr + l` for layered).
  /// Only "owned" (non-ghost) atlas nodes are stored, in the order they appear
  /// in the geometry's NodeColumns function space. Ranks load only their own
  /// nodes via a root-side read + scatter; root reassembles globally on flush.
  // -------------------------------------------------------------------------
  struct ThermoFrame {
    std::size_t ncat   = 0;
    std::size_t iceLev = 0;
    std::size_t snoLev = 0;
    std::size_t nOwnedNodes = 0;
    std::vector<double> Tsfcn;  // [node, k]
    std::vector<double> qice;   // [node, k, l] over iceLev
    std::vector<double> sice;   // [node, k, l] over iceLev
    std::vector<double> qsno;   // [node, k, l] over snoLev
    std::vector<double> apnd;   // [node, k]
    std::vector<double> hpnd;   // [node, k]
    std::vector<double> ipnd;   // [node, k]

    // Convenience accessors.
    double & at2(std::vector<double> & a,
                 std::size_t node, std::size_t k) const {
      return a[node * ncat + k];
    }
    double at2(const std::vector<double> & a,
               std::size_t node, std::size_t k) const {
      return a[node * ncat + k];
    }
    double & at3(std::vector<double> & a, std::size_t nlyr,
                 std::size_t node, std::size_t k, std::size_t l) const {
      return a[(node * ncat + k) * nlyr + l];
    }
    double at3(const std::vector<double> & a, std::size_t nlyr,
               std::size_t node, std::size_t k, std::size_t l) const {
      return a[(node * ncat + k) * nlyr + l];
    }
  };

  /// Read the thermo + pond variables from the input restart on root, scatter
  /// to all ranks. `ncat`, `iceLev`, `snoLev` must match the input file.
  ThermoFrame readThermo(std::size_t ncat,
                         std::size_t iceLev,
                         std::size_t snoLev) const;

  /// Gather the modified thermo + pond data to root and overwrite the
  /// corresponding variables in the output restart (which must already exist
  /// from a prior `write()` of the FieldSet).
  void flushThermo(const ThermoFrame & frame) const;

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

  /// Build the list of owned (gidx) for the local rank. Returns sorted
  /// vector of 1-based global indices in node-iteration order.
  std::vector<atlas::gidx_t> ownedGlobalIndices() const;
};

}  // namespace soca
