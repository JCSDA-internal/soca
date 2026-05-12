/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/array.h"
#include "atlas/field/FieldSet.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/ObjectHandle.h"
#include "eckit/config/Configuration.h"
#include "oops/util/Printable.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "soca/PostProcess/CiceRestartIO.h"

namespace soca {

class Geometry;
class State;

class PostProcessIce: public util::Printable {
 public:
  // ---------------------------------------------------------------------------
  /// @brief Parameters for the SST update on ice2noice transitions. When ice
  /// is removed by the analysis, optionally warm the surface ocean toward the
  /// freezing temperature so the post-DA SST is consistent with no ice.
  class SstUpdateParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(SstUpdateParameters, oops::Parameters)
   public:
    oops::Parameter<bool> adjustSST{"update SST",
      "Warm the surface ocean toward Tfrz on ice2noice transitions.",
      true, this};
    oops::Parameter<double> sstDiffMax{"max positive SST update",
      "Maximum update to sea surface temperature (K).", 1.0, this,
      {oops::minConstraint(0.0)}};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for re-binning the ice thickness distribution after
  /// rescale, so that each category's mean thickness is inside its bin.
  class ITDParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(ITDParameters, oops::Parameters)
   public:
    oops::Parameter<bool> rebin{"rebin",
      "If true, re-bin ice thickness categories so that hin = vicen/aicen "
      "is inside [hicat[k], hicat[k+1]] for every category.", false, this};
    oops::Parameter<std::vector<double>> hicat{"category bounds",
      "Lower edges of CICE thickness categories plus the upper edge of the "
      "thickest category. Length must equal ncat+1.",
      {0.0, 0.6445072, 1.391433, 2.470179, 4.567288, 9.333887}, this};
    oops::Parameter<double> dhiMin{"min thickness gap",
      "Minimum gap from category bounds enforced during re-binning (m).",
      0.01, this, {oops::minConstraint(0.0)}};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters controlling the snow-volume distribution. The cell-mean
  /// snow thickness is clamped to [hsnowMin, hsnowMax] and then distributed
  /// per category proportional to aicen (matches the Python reference
  /// scripts' insert_hsnow path).
  class SnowParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(SnowParameters, oops::Parameters)
   public:
    oops::Parameter<double> hsnowMin{"min snow thickness",
      "Floor on per-category snow thickness (m). Categories with vsnon/aicen "
      "below this are zeroed after rescale. Set to 0 to disable.",
      0.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> hsnowMax{"max snow thickness",
      "Cap on the cell-mean snow thickness (m). The analysis hsno is clamped "
      "to [min snow thickness / aice, this] before being distributed per cat. "
      "Matches the Python reference scripts. Set very large to disable.",
      5.0, this, {oops::minConstraint(0.0)}};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for the snow / ice freeboard enforcement step. After
  /// snow has been redistributed, each category should satisfy the hydrostatic
  /// balance rho_ice*hi + rho_snow*hs <= rho_ocean*(hi+hs). Snow is first
  /// redistributed across categories; if any cat is still flooded, ice volume
  /// is grown to lift the snow-ice interface back to sea level.
  class FreeboardParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(FreeboardParameters, oops::Parameters)
   public:
    oops::Parameter<bool> enforce{"enforce",
      "If true, run the freeboard enforcement pass after the ITD rebin.",
      false, this};
    oops::Parameter<double> rhoIce{"rho ice",
      "Sea ice density (kg/m^3)", 917.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> rhoSnow{"rho snow",
      "Snow density (kg/m^3)",   330.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> rhoOcean{"rho ocean",
      "Ocean density (kg/m^3)", 1025.0, this, {oops::minConstraint(0.0)}};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for the Stage C thermo / pond pass. Operates on the CICE
  /// restart's per-layer enthalpy, surface temperature, and melt-pond fields
  /// (read directly via CiceRestartIO::readThermo, since these variables are
  /// not part of the soca State / FieldSet pipeline).
  class ThermoParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(ThermoParameters, oops::Parameters)

   public:
    oops::Parameter<bool> updateSnowThermo{"update snow thermo",
      "Clip snow enthalpy into [snowEnthalpy(min Tsfc), snowEnthalpy(max Tsfc)] "
      "and back-derive Tsfcn from qsno on cats with snow. Also caps the surface "
      "ice layer enthalpy by iceEnthalpyBL99(Tsfcn, sice). On by default to "
      "match the Python reference scripts.",
      true, this};
    oops::Parameter<bool> resetPonds{"reset ponds",
      "Zero apnd and hpnd on per-cat slots where the snow distribution "
      "modified vsnon. Mirrors the Python reference, which only zeros ponds "
      "on snow-inserted cats. Never touches ipnd. On by default.",
      true, this};
    oops::Parameter<double> maxTsfc{"max Tsfc",
      "Upper bound on Tsfcn (deg C). Snow enthalpy is clipped accordingly.",
      -1.0, this};
    oops::Parameter<double> minTsfc{"min Tsfc",
      "Lower bound on Tsfcn (deg C). Snow enthalpy is clipped accordingly.",
      -100.0, this};
    oops::Parameter<bool> seedNewIce{"seed new ice",
      "Seed Tsfcn/qsno/sice/qice for cats that went from aicen=0 in background "
      "to aicen>0 after Stage A/B. Donor Tsfc from a global lat/lon nearest "
      "neighbor with any ice; sub-surface profile from CICE physics "
      "(BL99/BZ99). On by default to match the Python reference scripts.",
      true, this};
    oops::Parameter<int> seedSearchK{"seed search neighbors",
      "Number of nearest KDTree neighbors to scan for a donor with any ice "
      "before falling back to Tfrz seeding.",
      64, this, {oops::minConstraint(1)}};
  };

  // ---------------------------------------------------------------------------
  /// @brief Input/output paths for the CICE restart that this postprocess pass
  /// reads from and writes to. Both required.
  class CiceRestartParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(CiceRestartParameters, oops::Parameters)
   public:
    oops::RequiredParameter<std::string> input{"input",
      "Path to the input CICE restart NetCDF file.", this};
    oops::RequiredParameter<std::string> output{"output",
      "Path to the output CICE restart NetCDF file (will be overwritten).", this};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for adding soca increment to CICE restart files
  class Parameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Parameters, oops::Parameters)

   public:
    oops::RequiredParameter<int> ncat{"ncat", this, {oops::minConstraint(1)}};
    oops::RequiredParameter<int> ice_lev{"ice_lev", this, {oops::minConstraint(1)}};
    oops::RequiredParameter<int> sno_lev{"sno_lev", this, {oops::minConstraint(1)}};
    oops::Parameter<SstUpdateParameters> sstUpdate{"sst update",
      "SST adjustment on ice2noice transitions.", {}, this};
    oops::Parameter<ITDParameters> itd{"itd", "ITD re-bin options", {}, this};
    oops::Parameter<SnowParameters> snow{"snow",
      "Snow-volume distribution options.", {}, this};
    oops::Parameter<FreeboardParameters> freeboard{"freeboard",
      "Freeboard enforcement options.", {}, this};
    oops::Parameter<ThermoParameters> thermo{"thermo",
      "Stage C thermo / pond options.", {}, this};
    oops::RequiredParameter<CiceRestartParameters> ciceRestart{"cice restart",
      "Input/output paths for the CICE restart this pass reads and writes.",
      this};
    oops::Parameter<double> minAice{"min aice",
            "minimum allowable ice concentration", 0.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> minVice{"min vice",
            "minimum allowable ice volume", 0.00001, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> minAiceOutput{"min aice output",
            "Cells where the analysis aice is below this threshold are treated "
            "as no-ice in the output (all per-cat fields zeroed). Stops "
            "near-zero aice (and the rebin's resulting micro-vicen) from "
            "polluting the CICE restart. Set to 0 to disable.",
            1.0e-4, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> hitotMin{"min new ice thickness",
            "Minimum cell-mean ice thickness (m) used to set the rebin volume "
            "target when noice->ice transitions occur and no KDTree donor with "
            "ice can be found within seed-search-neighbors. Mirrors the Python "
            "reference scripts' hitot_min.",
            0.1, this, {oops::minConstraint(0.0)}};
  };

  const std::string classname() {return "soca::PostProcessIce";}

  PostProcessIce(const Geometry &, const eckit::Configuration &);

  void postProcess(State & pproc, const State & restart,
                   const State & analysis) const;

  /// Stage C thermo / pond pass. Mutates `frame` in place using the per-cat
  /// aicen / vsnon fields from `fset`. Owned-node ordering must match the
  /// CiceRestartIO traversal (ghost == 0 && global_index > 0). No-op when
  /// neither `update snow thermo` nor `reset ponds` is enabled.
  ///
  /// `snowTouched` is sized `nOwnedNodes * ncat`, indexed `[on * ncat + k]`,
  /// non-zero where the snow distribution step modified vsnon on that
  /// (owned-node, cat) slot. Used to gate the apnd/hpnd reset; ipnd is never
  /// reset. May be empty when `reset ponds` is off.
  void applyThermoStage(CiceRestartIO::ThermoFrame & frame,
                        const atlas::FieldSet & fset,
                        const std::vector<std::uint8_t> & snowTouched) const;

  /// @brief Per-cell donor record assembled by the sparse halo exchange. Holds
  /// the donor cell's per-category ice and thermo values so consumers
  /// (donorMeanIce in the NOICE2ICE branch, seedNewIce in Stage C) can read
  /// from them directly.
  ///
  /// All four arrays are flat, sized in terms of `ncat`, `iceLev`, `snoLev`
  /// from `params_`. Layout:
  ///   aicen, vicen, vsnon, Tsfcn: indexed [k]
  ///   qice, sice:                 indexed [k * iceLev + l]
  ///   qsno:                       indexed [k * snoLev + l]
  struct CatRecord {
    std::vector<double> aicen;       // ncat
    std::vector<double> vicen;       // ncat
    std::vector<double> vsnon;       // ncat
    std::vector<double> Tsfcn;       // ncat
    std::vector<double> qice;        // ncat * iceLev
    std::vector<double> sice;        // ncat * iceLev
    std::vector<double> qsno;        // ncat * snoLev
    double mask;
  };

 private:
  const Geometry & geom_;
  Parameters params_;
  size_t ncat_;
  // Global lat/lon KDTree. Payload is the 1-based atlas global_index (gidx).
  // Built at construction by all-gathering owned-cell (lon, lat, gidx) triples
  // across all ranks and building one tree on each rank.
  atlas::util::IndexKDTree kdTree_;
  // Owner rank for each gidx in the global cell set, derived from the per-rank
  // counts during the kdTree allGatherv. Used to route donor-data requests in
  // the sparse halo exchange.
  std::unordered_map<std::int64_t, int> gidxToOwnerRank_;
  // Local-rank only: owned-node atlas jnode for each owned gidx. Used when
  // packing reply records for donor-data requests received from other ranks.
  std::unordered_map<std::int64_t, std::size_t> gidxToLocalJnode_;

  atlas::array::ArrayView<double, 2> lonlat_;
  atlas::array::ArrayView<double, 2> mask_;

  /// @brief Mask of (ownedNode, k) cells that transitioned from
  /// `bg_aicen[k] == 0` to `new_aicen[k] > 0` after Stages A/B. Indexed
  /// `ownedNode * ncat + k`.
  struct NewIceMask {
    std::size_t ncat = 0;
    std::size_t nOwnedNodes = 0;
    std::vector<std::uint8_t> data;
    bool at(std::size_t ownedNode, std::size_t k) const {
      return data[ownedNode * ncat + k] != 0;
    }
  };

  /// Phase A+B+C of the sparse halo exchange. For every owned cell, run
  /// `kdTree_.closestPoints(target, K)` to identify the global gidx of K
  /// nearest neighbors; deduplicate into a per-rank "wanted" set; sparse
  /// allToAllv to request donor data; pack reply records from local FieldSet
  /// views and the local thermo frame; sparse allToAllv to receive replies.
  /// Returns a map keyed by gidx with the donor's full CatRecord.
  std::unordered_map<std::int64_t, CatRecord> gatherDonorHalo(
      std::size_t K,
      const std::vector<atlas::array::ArrayView<double, 2>> & bg_aice_cat,
      const std::vector<atlas::array::ArrayView<double, 2>> & bg_vice_cat,
      const std::vector<atlas::array::ArrayView<double, 2>> & bg_vsno_cat,
      const CiceRestartIO::ThermoFrame & frame,
      const std::vector<std::int64_t> & ownedNodeOf,
      std::size_t ice_lev,
      std::size_t sno_lev) const;

  /// Area-weighted mean Tsfc and mean ice thickness across the K nearest
  /// KDTree neighbors that carry any ice. Used by the Stage A case-dispatch
  /// for noice→ice cells to set the rebin's volume target before the per-
  /// category thermo seeding runs in Stage C.
  ///
  /// Returns true if at least one donor with ice was found. Outputs are valid
  /// only when the return is true; otherwise both are left untouched.
  bool donorMeanIce(std::size_t jnode,
                    std::size_t K,
                    const std::unordered_map<std::int64_t, CatRecord> &
                        donorCache,
                    double & Tsfc_mean,
                    double & hice_mean) const;

  /// Stage C noice→ice seeding. For each (ownedNode, k) flagged in `mask`,
  /// pick the global lat/lon nearest cell with any ice from `donorCache` and
  /// seed Tsfcn from its area-weighted mean Tsfc; synthesize qsno/qice/sice
  /// from CICE physics (snowEnthalpy, iceEnthalpyBL99, siceLayerCice4).
  /// Falls back to Tfrz physics seed when no donor is found within K.
  /// Returns the count of fall-back cells for logging.
  std::size_t seedNewIce(CiceRestartIO::ThermoFrame & frame,
                         const NewIceMask & mask,
                         const std::unordered_map<std::int64_t, CatRecord> & donorCache,
                         const std::vector<std::int64_t> & ownedNodeOf,
                         std::size_t ice_lev,
                         std::size_t sno_lev) const;

  // Helpers to compute totals and mean thicknesses from category fields
  double totalAice(const std::vector<atlas::array::ArrayView<double, 2>> & aiceCat,
                   size_t jnode) const;
  double meanHice(const std::vector<atlas::array::ArrayView<double, 2>> & viceCat,
                  double aice, size_t jnode) const;
  double meanHsno(const std::vector<atlas::array::ArrayView<double, 2>> & vsnoCat,
                  double aice, size_t jnode) const;
  void print(std::ostream &) const override;
};

}  // namespace soca
