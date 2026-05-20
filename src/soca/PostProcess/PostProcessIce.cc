/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/PostProcess/PostProcessIce.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_set>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/util/Point.h"

#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/PostProcess/IcePhysics.h"
#include "soca/State/State.h"

namespace soca {

// -----------------------------------------------------------------------------
// Helpers to gather atlas array views for the per-category and
// per-(category,layer) CICE thermo/pond fields. Field names follow the
// convention declared in fields_metadata.yml:
//   per-cat:        <prefix><c><suffix>           e.g. sea_ice_category3_pond_depth
//   per-cat-layer:  <prefix><c><mid><l><suffix>   e.g. sea_ice_category3_layer5_enthalpy
// Category and layer indices are 1-based in field names.

namespace {

std::vector<atlas::array::ArrayView<double, 2>> ppiCatViews(
    const atlas::FieldSet & fset, std::size_t ncat,
    const std::string & prefix, const std::string & suffix) {
  std::vector<atlas::array::ArrayView<double, 2>> views;
  views.reserve(ncat);
  for (std::size_t k = 0; k < ncat; ++k) {
    const std::string name = prefix + std::to_string(k + 1) + suffix;
    views.push_back(atlas::array::make_view<double, 2>(fset.field(name)));
  }
  return views;
}

std::vector<std::vector<atlas::array::ArrayView<double, 2>>> ppiCatLevViews(
    const atlas::FieldSet & fset, std::size_t ncat, std::size_t nlev,
    const std::string & prefix, const std::string & mid,
    const std::string & suffix) {
  std::vector<std::vector<atlas::array::ArrayView<double, 2>>> views(ncat);
  for (std::size_t k = 0; k < ncat; ++k) {
    views[k].reserve(nlev);
    for (std::size_t l = 0; l < nlev; ++l) {
      const std::string name = prefix + std::to_string(k + 1) + mid
                             + std::to_string(l + 1) + suffix;
      views[k].push_back(atlas::array::make_view<double, 2>(fset.field(name)));
    }
  }
  return views;
}

}  // namespace

// -----------------------------------------------------------------------------

PostProcessIce::PostProcessIce(const Geometry & geom,
                               const eckit::Configuration & conf)
  : geom_(geom),
    lonlat_(atlas::array::make_view<double, 2>(geom.functionSpace().lonlat())),
    mask_(atlas::array::make_view<double, 2>(geom.fields().field("interp_mask"))) {
  params_.deserialize(conf);
  ncat_ = params_.ncat.value();

  // Build a global lat/lon KDTree across all ranks. Payload is the 1-based
  // atlas global_index (gidx). Pattern adapted from oops::GeometryData::setGlobalTree.
  const auto & fs   = geom_.functionSpace();
  const auto ghost  = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindex = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());
  std::vector<double> ownedLon, ownedLat;
  std::vector<atlas::gidx_t> ownedGidx;
  std::vector<std::size_t> ownedJnode;
  ownedLon.reserve(ghost.size());
  ownedLat.reserve(ghost.size());
  ownedGidx.reserve(ghost.size());
  ownedJnode.reserve(ghost.size());
  for (std::size_t jnode = 0; jnode < ghost.size(); ++jnode) {
    if (ghost(jnode)) continue;
    if (gindex(jnode) <= 0) continue;
    ownedLon.push_back(lonlat_(jnode, 0));
    ownedLat.push_back(lonlat_(jnode, 1));
    ownedGidx.push_back(gindex(jnode));
    ownedJnode.push_back(jnode);
  }
  for (std::size_t i = 0; i < ownedGidx.size(); ++i) {
    gidxToLocalJnode_.emplace(static_cast<std::int64_t>(ownedGidx[i]),
                              ownedJnode[i]);
  }

  const auto & comm = geom_.getComm();
  const std::size_t nr = comm.size();
  std::vector<std::size_t> sizes(nr, 0);
  comm.allGather(ownedLon.size(), sizes.begin(), sizes.end());

  std::vector<double> globalLon, globalLat;
  std::vector<atlas::gidx_t> globalGidx;
  oops::mpi::allGatherv(comm, ownedLon, globalLon);
  oops::mpi::allGatherv(comm, ownedLat, globalLat);
  oops::mpi::allGatherv(comm, ownedGidx, globalGidx);

  ASSERT(globalLon.size() == globalLat.size());
  ASSERT(globalLon.size() == globalGidx.size());

  // Build owner-rank map from per-rank counts (concatenation order matches
  // allGatherv contract).
  std::size_t cursor = 0;
  for (std::size_t r = 0; r < nr; ++r) {
    for (std::size_t i = 0; i < sizes[r]; ++i) {
      gidxToOwnerRank_.emplace(static_cast<std::int64_t>(globalGidx[cursor]),
                               static_cast<int>(r));
      ++cursor;
    }
  }
  ASSERT(cursor == globalGidx.size());

  // Build the KDTree with gidx as payload.
  std::vector<atlas::PointLonLat> nodes;
  std::vector<atlas::gidx_t> payloads;
  nodes.reserve(globalLon.size());
  payloads.reserve(globalLon.size());
  for (std::size_t i = 0; i < globalLon.size(); ++i) {
    nodes.emplace_back(globalLon[i], globalLat[i]);
    payloads.push_back(globalGidx[i]);
  }
  kdTree_.build(nodes, payloads);
}

// -----------------------------------------------------------------------------

void PostProcessIce::postProcess(State & pproc,
                                 const State & restart,
                                 const State & analysis) const {
  // Access function space and sizes
  const auto & fs = geom_.functionSpace();
  // Restart background fields (read-only)
  const atlas::FieldSet & bgfields = restart.fieldSet();
  // Analysis fields (read-only)
  atlas::FieldSet anfields = analysis.fieldSet();
  pproc = restart;  // start from background state
  atlas::FieldSet & newfields = pproc.fieldSet();
  oops::Log::info() << " background restart: " << restart << std::endl;
  oops::Log::info() << " analysis: " << analysis << std::endl;
  oops::Log::info() << " pproc before " << pproc << std::endl;

  const size_t ice_lev = params_.ice_lev.value();
  const size_t sno_lev = params_.sno_lev.value();
  const size_t field_size = bgfields.field(0).shape(0);

  // Restart background fields
  std::vector<atlas::array::ArrayView<double, 2>> bg_aice_cat, bg_vice_cat, bg_vsno_cat;
  std::vector<atlas::array::ArrayView<double, 2>> new_aice_cat, new_vice_cat, new_vsno_cat;
  bg_aice_cat.reserve(ncat_);
  bg_vice_cat.reserve(ncat_);
  bg_vsno_cat.reserve(ncat_);
  new_aice_cat.reserve(ncat_);
  new_vice_cat.reserve(ncat_);
  new_vsno_cat.reserve(ncat_);
  std::string varname;
  for (size_t icat = 0; icat < ncat_; ++icat) {
    varname = "sea_ice_category" + std::to_string(icat+1) + "_area_fraction";
    bg_aice_cat.push_back(atlas::array::make_view<double, 2>(bgfields.field(varname)));
    new_aice_cat.push_back(atlas::array::make_view<double, 2>(newfields.field(varname)));
    varname = "sea_ice_category" + std::to_string(icat+1) + "_volume";
    bg_vice_cat.push_back(atlas::array::make_view<double, 2>(bgfields.field(varname)));
    new_vice_cat.push_back(atlas::array::make_view<double, 2>(newfields.field(varname)));
    varname = "sea_ice_snow_category" + std::to_string(icat+1) + "_volume";
    bg_vsno_cat.push_back(atlas::array::make_view<double, 2>(bgfields.field(varname)));
    new_vsno_cat.push_back(atlas::array::make_view<double, 2>(newfields.field(varname)));
  }
  // Writable view onto pproc's ocean temperature so the SST-adjust step on
  // ice2noice cells can update it (Fortran soca2cice behavior).
  auto new_tocn = atlas::array::make_view<double, 2>(
      newfields.field("sea_water_potential_temperature"));
  // Create field and view for total ice concentration in restarts
  atlas::Field bg_aice_field = fs.createField<double>(
    atlas::option::name("sea_ice_area_fraction") | atlas::option::levels(1));
  atlas::Field new_aice_field = fs.createField<double>(
    atlas::option::name("sea_ice_area_fraction") | atlas::option::levels(1));
  atlas::Field new_hice_field = fs.createField<double>(
    atlas::option::name("sea_ice_thickness") | atlas::option::levels(1));
  atlas::Field new_hsno_field = fs.createField<double>(
    atlas::option::name("sea_ice_snow_thickness") | atlas::option::levels(1));
  auto bg_aice = atlas::array::make_view<double, 2>(bg_aice_field);
  auto new_aice = atlas::array::make_view<double, 2>(new_aice_field);
  auto new_hice = atlas::array::make_view<double, 2>(new_hice_field);
  auto new_hsno = atlas::array::make_view<double, 2>(new_hsno_field);
  newfields.add(new_aice_field);
  newfields.add(new_hice_field);
  newfields.add(new_hsno_field);
  // Compute total ice concentration in restart
  for (size_t jnode = 0; jnode < field_size; ++jnode) {
    bg_aice(jnode, 0) = totalAice(bg_aice_cat, jnode);
    new_aice(jnode, 0) = bg_aice(jnode, 0);
    new_hice(jnode, 0) = meanHice(new_vice_cat, new_aice(jnode, 0), jnode);
    new_hsno(jnode, 0) = meanHsno(new_vsno_cat, new_aice(jnode, 0), jnode);
  }

  // Analysis fields. The user lists which ice-related variables are actually
  // produced by the DA in `params_.analysisVars`. We resolve each of the three
  // logical inputs (aice, hice, hsno) from that list:
  //   - aice: required; must be `sea_ice_area_fraction`.
  //   - ice volume: either `sea_ice_thickness` (per-ice-area, m) or
  //                 `sea_ice_volume` (per-cell-area, m). Missing -> bg fallback.
  //   - snow volume: either `sea_ice_snow_thickness` or `sea_ice_snow_volume`.
  //                  Missing -> bg fallback (the per-cell loop already does
  //                  this when its per-jnode target is non-positive).
  // Internally we always work in per-ice-area thickness; volume forms are
  // divided by aice on read.
  const oops::Variables & anVars = params_.analysisVars.value();
  if (!anVars.has("sea_ice_area_fraction")) {
    throw eckit::UserError(
      "PostProcessIce: `analysis variables` must include sea_ice_area_fraction",
      Here());
  }
  const bool has_hice_th  = anVars.has("sea_ice_thickness");
  const bool has_hice_vol = anVars.has("sea_ice_volume");
  const bool has_hsno_th  = anVars.has("sea_ice_snow_thickness");
  const bool has_hsno_vol = anVars.has("sea_ice_snow_volume");
  if (has_hice_th && has_hice_vol) {
    throw eckit::UserError(
      "PostProcessIce: `analysis variables` lists both sea_ice_thickness and "
      "sea_ice_volume; pick one", Here());
  }
  if (has_hsno_th && has_hsno_vol) {
    throw eckit::UserError(
      "PostProcessIce: `analysis variables` lists both sea_ice_snow_thickness "
      "and sea_ice_snow_volume; pick one", Here());
  }

  auto a_aice = atlas::array::make_view<double, 2>(
      anfields.field("sea_ice_area_fraction"));

  // Per-jnode per-ice-area thickness targets, materialized once. Negative-as-
  // sentinel "no analysis here, fall back to bg" propagates through the
  // existing per-cell loop logic.
  std::vector<double> a_hice_vec(field_size, -1.0);
  std::vector<double> a_hsno_vec(field_size, -1.0);
  if (has_hice_th) {
    auto v = atlas::array::make_view<double, 2>(
        anfields.field("sea_ice_thickness"));
    for (size_t jn = 0; jn < field_size; ++jn) a_hice_vec[jn] = v(jn, 0);
  } else if (has_hice_vol) {
    auto v = atlas::array::make_view<double, 2>(
        anfields.field("sea_ice_volume"));
    for (size_t jn = 0; jn < field_size; ++jn) {
      const double ai = a_aice(jn, 0);
      a_hice_vec[jn] = (ai > icephysics::Constants::puny) ? v(jn, 0) / ai : 0.0;
    }
  } else {
    // Bg fallback: per-ice-area thickness from background per-cat volumes.
    for (size_t jn = 0; jn < field_size; ++jn) {
      double v = 0.0;
      for (size_t k = 0; k < ncat_; ++k) v += bg_vice_cat[k](jn, 0);
      a_hice_vec[jn] = (bg_aice(jn, 0) > icephysics::Constants::puny)
                       ? v / bg_aice(jn, 0) : 0.0;
    }
  }
  if (has_hsno_th) {
    auto v = atlas::array::make_view<double, 2>(
        anfields.field("sea_ice_snow_thickness"));
    for (size_t jn = 0; jn < field_size; ++jn) a_hsno_vec[jn] = v(jn, 0);
  } else if (has_hsno_vol) {
    auto v = atlas::array::make_view<double, 2>(
        anfields.field("sea_ice_snow_volume"));
    for (size_t jn = 0; jn < field_size; ++jn) {
      const double ai = a_aice(jn, 0);
      a_hsno_vec[jn] = (ai > icephysics::Constants::puny) ? v(jn, 0) / ai : 0.0;
    }
  }
  // hsno bg-fallback is left to the existing per-cell snow block: when
  // a_hsno_vec[jnode] is non-positive (the default sentinel) it recovers
  // hsn_target from the bg snow cell-mean. So nothing to do here.

  // parameters
  const auto & sstUpdateParams = params_.sstUpdate.value();
  const bool   adjust_sst = sstUpdateParams.adjustSST.value();
  const double max_dsst   = sstUpdateParams.sstDiffMax.value();
  const double min_aice   = params_.minAice.value();
  const double min_vice = params_.minVice.value();
  const auto & itd = params_.itd.value();
  const bool do_rebin = itd.rebin.value();
  const std::vector<double> hicat = itd.hicat.value();
  const double dhi_min = itd.dhiMin.value();
  if (do_rebin && hicat.size() != ncat_ + 1) {
    throw eckit::UserError(
      "PostProcessIce: itd.category bounds must have length ncat+1", Here());
  }
  // Per-cell scratch buffers for the ITD rebin.
  std::vector<double> rebin_aicen(ncat_, 0.0);
  std::vector<double> rebin_vicen(ncat_, 0.0);
  size_t rebin_failures = 0;

  const auto & snowParams = params_.snow.value();
  const double hsnow_min = snowParams.hsnowMin.value();

  const auto & fbParams = params_.freeboard.value();
  const bool do_freeboard = fbParams.enforce.value();
  const double rho_ice   = fbParams.rhoIce.value();
  const double rho_snow  = fbParams.rhoSnow.value();
  const double rho_ocean = fbParams.rhoOcean.value();
  // Per-cell scratch buffers for freeboard.
  std::vector<double> fb_aicen(ncat_, 0.0);
  std::vector<double> fb_vicen(ncat_, 0.0);
  std::vector<double> fb_vsnon(ncat_, 0.0);
  size_t freeboard_failures = 0;

  // ---------------------------------------------------------------------------
  // CICE thermo/pond background views. These fields are carried as atlas Fields
  // of the `restart` State (declared in fields_metadata.yml, listed in the
  // State's `state variables`). `pproc = restart` already copied them into
  // `pproc.fieldSet()`; the background copies below are read for the donor
  // halo gather, while the per-cell loop and Stage C mutate `pproc.fieldSet()`.
  const auto bg_tsfc_cat     = ppiCatViews(restart.fieldSet(), ncat_,
                                           "sea_ice_category", "_surface_temperature");
  const auto bg_qice_cat_lev = ppiCatLevViews(restart.fieldSet(), ncat_, ice_lev,
                                              "sea_ice_category", "_layer", "_enthalpy");
  const auto bg_sice_cat_lev = ppiCatLevViews(restart.fieldSet(), ncat_, ice_lev,
                                              "sea_ice_category", "_layer", "_ice_salinity");
  const auto bg_qsno_cat_lev = ppiCatLevViews(restart.fieldSet(), ncat_, sno_lev,
                                              "sea_ice_snow_category", "_layer", "_enthalpy");

  // Sparse halo exchange for donor data (Phases A+B+C). K is the seed-new-ice
  // search radius; both donorMeanIce (in the NOICE2ICE branch) and seedNewIce
  // (Stage C) read from the same donorCache.
  const std::size_t seedK = static_cast<std::size_t>(
      params_.thermo.value().seedSearchK.value());
  auto donorCache = gatherDonorHalo(seedK, bg_aice_cat, bg_vice_cat, bg_vsno_cat,
                                    bg_tsfc_cat, bg_qice_cat_lev, bg_sice_cat_lev,
                                    bg_qsno_cat_lev, ice_lev, sno_lev);

  // Additional knobs picked up from the top-level parameters. seedK was
  // already computed above for the donor halo gather; reuse it here.
  const double hsnow_max        = snowParams.hsnowMax.value();
  const double hitot_min        = params_.hitotMin.value();
  const double min_aice_output  = params_.minAiceOutput.value();

  // Per-cell scratch for the rebin call.
  int      rebin_aicen_mutated   = 0;
  double   rebin_max_delta_aicen = 0.0;
  std::size_t rebin_aicen_mutations  = 0;
  double   rebin_max_delta_aicen_all = 0.0;
  size_t rebin_visited = 0;

  // ----------------------------------------------------------------------
  // Stage A: per-cell case dispatch (port of the Python reference scripts).
  // For each owned cell:
  //   1. clamp analysis bounds;
  //   2. classify (LAND/ICE2NOICE/NOICE2ICE/ICE2ICE);
  //   3. set seed aicen/vicen and vtot_target accordingly;
  //   4. ITD rebin (when enabled) drives aicen/vicen to (a_aice, vtot_target);
  //   5. distribute snow by aicen-weight from clamped cell-mean hsno;
  //   6. freeboard enforcement (when enabled);
  //   7. min-vice cleanup; recompute aggregates.
  for (size_t jnode = 0; jnode < field_size; ++jnode) {
    // Clamp bounds on analysis first.
    if (a_aice(jnode, 0) < min_aice) a_aice(jnode, 0) = 0.0;
    if (a_aice(jnode, 0) > 1.0)      a_aice(jnode, 0) = 1.0;
    if (a_hice_vec[jnode] < 0.0)     a_hice_vec[jnode] = 0.0;
    if (a_hsno_vec[jnode] < 0.0)     a_hsno_vec[jnode] = 0.0;
    // Drop output ice on cells where analysis aice is below the noise floor.
    // The rebin can otherwise produce micro-vicen layouts that show up as
    // thick-ice-on-edge artifacts after the min_vice cleanup.
    if (a_aice(jnode, 0) > 0.0 && a_aice(jnode, 0) < min_aice_output) {
      a_aice(jnode, 0) = 0.0;
    }

    const double ai_bg = bg_aice(jnode, 0);
    const double ai_an = a_aice(jnode, 0);
    const bool   isLand = (mask_(jnode, 0) == 0.0);

    // -------------------- Case 1: LAND or ICE2NOICE -----------------------
    if (isLand || ai_an <= icephysics::Constants::puny) {
      const bool wasIce = !isLand && (ai_bg > icephysics::Constants::puny);
      for (size_t icat = 0; icat < ncat_; ++icat) {
        new_aice_cat[icat](jnode, 0) = 0.0;
        new_vice_cat[icat](jnode, 0) = 0.0;
        new_vsno_cat[icat](jnode, 0) = 0.0;
      }
      // SST adjustment on ice2noice: warm the surface ocean by up to max_dsst
      // toward freezing (mirrors the Fortran soca2cice path).
      if (wasIce && adjust_sst) {
        const double Tfrz_loc = icephysics::Constants::Tfrz;
        const double dT = std::min(max_dsst,
                                   std::max(0.0, Tfrz_loc - new_tocn(jnode, 0)));
        new_tocn(jnode, 0) += dT;
      }
      new_aice(jnode, 0) = 0.0;
      new_hice(jnode, 0) = 0.0;
      new_hsno(jnode, 0) = 0.0;
      continue;
    }

    // -------------------- Build seed aicen/vicen --------------------------
    // Two cases: NOICE2ICE and ICE2ICE.
    double vtot_target = 0.0;
    bool   noice2ice   = (ai_bg <= icephysics::Constants::puny);

    if (noice2ice) {
      // Placeholder all-in-cat-0; the rebin will spread it. Vtot_target source
      // order:
      //   1. analysis a_hice (most informative — use it whenever positive);
      //   2. KDTree donor mean hice;
      //   3. hitot_min · ai_an fallback.
      // The donor-mean Tsfc is still useful for the Stage C thermo seed (it
      // sets initial Tsfcn/qsno/qice/sice physics consistently), so the
      // donorMeanIce call stays in the seedNewIce pass.
      double hice_seed = a_hice_vec[jnode];
      if (hice_seed <= 0.0) {
        double Tsfc_donor, hice_donor;
        const bool foundDonor = donorMeanIce(jnode, seedK, donorCache,
                                             Tsfc_donor, hice_donor);
        hice_seed = (foundDonor && hice_donor > hitot_min)
                        ? hice_donor : hitot_min;
      }
      vtot_target = hice_seed * ai_an;
      for (std::size_t k = 0; k < ncat_; ++k) {
        new_aice_cat[k](jnode, 0) = (k == 0) ? ai_an : 0.0;
        new_vice_cat[k](jnode, 0) = (k == 0) ? vtot_target : 0.0;
        new_vsno_cat[k](jnode, 0) = 0.0;
      }
    } else {
      // ICE2ICE: alpha-rescale per-cat from background. The rebin (if
      // enabled) then corrects toward analysis hice·a_aice. Trust the
      // analysis hice unconditionally (user decision 2026-05-11).
      const double alpha = ai_an / ai_bg;
      for (std::size_t k = 0; k < ncat_; ++k) {
        new_aice_cat[k](jnode, 0) = alpha * bg_aice_cat[k](jnode, 0);
        new_vice_cat[k](jnode, 0) = alpha * bg_vice_cat[k](jnode, 0);
        // vsnon is handled in the snow block below; clear here so the freeboard
        // pass doesn't see stale background snow.
        new_vsno_cat[k](jnode, 0) = 0.0;
      }
      vtot_target = a_hice_vec[jnode] * ai_an;
    }

    // -------------------- Stage B: ITD rebin ------------------------------
    if (do_rebin) {
      for (size_t icat = 0; icat < ncat_; ++icat) {
        rebin_aicen[icat] = new_aice_cat[icat](jnode, 0);
        rebin_vicen[icat] = new_vice_cat[icat](jnode, 0);
      }
      rebin_aicen_mutated   = 0;
      rebin_max_delta_aicen = 0.0;
      const bool ok = icephysics::adjustThicknessCategories(
          rebin_aicen, rebin_vicen, ai_an, vtot_target, hicat, dhi_min,
          /*ainMin=*/1.0e-8, /*aTol=*/1.0e-8, /*vTol=*/1.0e-6,
          &rebin_aicen_mutated, &rebin_max_delta_aicen);
      if (ok) {
        for (size_t icat = 0; icat < ncat_; ++icat) {
          new_aice_cat[icat](jnode, 0) = rebin_aicen[icat];
          new_vice_cat[icat](jnode, 0) = rebin_vicen[icat];
        }
        ++rebin_visited;
        if (rebin_aicen_mutated) ++rebin_aicen_mutations;
        rebin_max_delta_aicen_all =
            std::max(rebin_max_delta_aicen_all, rebin_max_delta_aicen);
      } else {
        ++rebin_failures;
      }
    }

    // -------------------- Stage C: snow distribution ----------------------
    // Recompute ai_now (may differ slightly from a_aice after rebin).
    double ai_now = 0.0;
    for (size_t icat = 0; icat < ncat_; ++icat) {
      ai_now += new_aice_cat[icat](jnode, 0);
    }
    if (ai_now > icephysics::Constants::puny) {
      // Target cell-mean snow thickness:
      //   * use analysis hsno when available;
      //   * fall back to background cell-mean.
      double hsn_target = a_hsno_vec[jnode];
      if (hsn_target <= 0.0) {
        // Recover from background.
        double bg_vsno_sum = 0.0;
        for (size_t icat = 0; icat < ncat_; ++icat) {
          bg_vsno_sum += bg_vsno_cat[icat](jnode, 0);
        }
        hsn_target = (bg_aice(jnode, 0) > icephysics::Constants::puny)
                         ? bg_vsno_sum / bg_aice(jnode, 0)
                         : 0.0;
      }
      // Clamp cell-mean to [hsnow_min, hsnow_max]. The min is treated as a
      // floor only when there's enough analysis to justify it; zero passes
      // through unchanged so noice/ice2noice still drop snow correctly.
      if (hsn_target > 0.0) {
        hsn_target = std::min(std::max(hsn_target, hsnow_min), hsnow_max);
      }
      const double vsno_total = hsn_target * ai_now;
      for (size_t icat = 0; icat < ncat_; ++icat) {
        new_vsno_cat[icat](jnode, 0) =
            vsno_total * (new_aice_cat[icat](jnode, 0) / ai_now);
      }
    }

    // -------------------- Stage C+: freeboard -----------------------------
    if (do_freeboard) {
      bool anyIce = false;
      for (size_t icat = 0; icat < ncat_; ++icat) {
        fb_aicen[icat] = new_aice_cat[icat](jnode, 0);
        fb_vicen[icat] = new_vice_cat[icat](jnode, 0);
        fb_vsnon[icat] = new_vsno_cat[icat](jnode, 0);
        if (fb_aicen[icat] > 0.0) anyIce = true;
      }
      if (anyIce) {
        const bool ok = icephysics::enforceFreeboard(
            fb_aicen, fb_vicen, fb_vsnon, rho_ice, rho_snow, rho_ocean);
        if (ok) {
          for (size_t icat = 0; icat < ncat_; ++icat) {
            new_aice_cat[icat](jnode, 0) = fb_aicen[icat];
            new_vice_cat[icat](jnode, 0) = fb_vicen[icat];
            new_vsno_cat[icat](jnode, 0) = fb_vsnon[icat];
          }
        } else {
          ++freeboard_failures;
        }
      }
    }

    // -------------------- min-vice cleanup (mass-conserving) --------------
    // Identify cats with vicen < min_vice. Sum their (aicen, vicen, vsnon)
    // and redistribute that mass into surviving cats proportionally to the
    // survivors' aicen, preserving Σaicen/Σvicen/Σvsnon. Without this, the
    // rebin's marginal-aice solutions (e.g. tiny mass in cat 0 + bulk in a
    // thicker cat) get clipped to one cat at its upper bin edge, producing
    // thick-ice-on-edge artifacts.
    {
      double dropA = 0.0, dropV = 0.0, dropS = 0.0;
      double survA = 0.0;
      for (size_t icat = 0; icat < ncat_; ++icat) {
        if (new_aice_cat[icat](jnode, 0) > 0.0
            && new_vice_cat[icat](jnode, 0) < min_vice) {
          dropA += new_aice_cat[icat](jnode, 0);
          dropV += new_vice_cat[icat](jnode, 0);
          dropS += new_vsno_cat[icat](jnode, 0);
          new_aice_cat[icat](jnode, 0) = 0.0;
          new_vice_cat[icat](jnode, 0) = 0.0;
          new_vsno_cat[icat](jnode, 0) = 0.0;
        } else {
          survA += new_aice_cat[icat](jnode, 0);
        }
      }
      if (dropA > 0.0 && survA > 0.0) {
        for (size_t icat = 0; icat < ncat_; ++icat) {
          if (new_aice_cat[icat](jnode, 0) > 0.0) {
            const double w = new_aice_cat[icat](jnode, 0) / survA;
            new_aice_cat[icat](jnode, 0) += dropA * w;
            new_vice_cat[icat](jnode, 0) += dropV * w;
            new_vsno_cat[icat](jnode, 0) += dropS * w;
          }
        }
      }
      // If no survivors, the dropped mass is gone — this happens only when
      // every cat had vicen < min_vice, i.e. the cell-mean vicen was below
      // the floor anyway. Acceptable.
    }
    // -------------------- Aggregate diagnostics ---------------------------
    new_aice(jnode, 0) = totalAice(new_aice_cat, jnode);
    new_hice(jnode, 0) = meanHice(new_vice_cat, new_aice(jnode, 0), jnode);
    new_hsno(jnode, 0) = meanHsno(new_vsno_cat, new_aice(jnode, 0), jnode);
  }
  if (do_rebin) {
    oops::Log::info() << "PostProcessIce: ITD rebin visited "
                      << rebin_visited << " cells; aicen-mutated "
                      << rebin_aicen_mutations << " (max |Δaicen|="
                      << rebin_max_delta_aicen_all << ")." << std::endl;
    if (rebin_failures > 0) {
      oops::Log::warning() << "PostProcessIce: ITD rebin failed at "
                           << rebin_failures << " cells (target outside the "
                           << "feasible envelope); left untouched." << std::endl;
    }
  }
  if (do_freeboard && freeboard_failures > 0) {
    oops::Log::warning() << "PostProcessIce: freeboard enforcement failed at "
                         << freeboard_failures << " cells; left untouched."
                         << std::endl;
  }
  oops::Log::info() << " after pp restart: " << restart << std::endl;

  // ---------------------------------------------------------------------------
  // Stage C: build NewIceMask, run applyThermoStage, optional seedNewIce.
  // With the shuffle path removed, every (jnode, k) cell that transitioned
  // bg_aicen==0 -> new_aicen>0 needs thermo seeding from a KDTree donor.
  // The masks are sized field_size*ncat and indexed by jnode (ghost cells
  // are simply left zero).
  NewIceMask newIce;
  newIce.ncat = ncat_;
  newIce.nNodes = field_size;
  newIce.data.assign(field_size * ncat_, 0);
  for (std::size_t jnode = 0; jnode < field_size; ++jnode) {
    for (std::size_t k = 0; k < ncat_; ++k) {
      const bool wasZero = (bg_aice_cat[k](jnode, 0) == 0.0);
      const bool nowPos  = (new_aice_cat[k](jnode, 0) > 0.0);
      if (wasZero && nowPos) {
        newIce.data[jnode * ncat_ + k] = 1;
      }
    }
  }

  const auto & thermoParams = params_.thermo.value();

  // Build snowTouched mask: non-zero on (jnode, k) slots where the snow
  // distribution modified vsnon vs the background restart. Drives the
  // apnd/hpnd reset in applyThermoStage; mirrors dd's policy of resetting
  // ponds only where snow was actually inserted.
  std::vector<std::uint8_t> snowTouched;
  snowTouched.assign(field_size * ncat_, 0);
  const double snowTouchedTol = 1.0e-12;
  for (std::size_t jnode = 0; jnode < field_size; ++jnode) {
    for (std::size_t k = 0; k < ncat_; ++k) {
      const double dv = std::fabs(new_vsno_cat[k](jnode, 0)
                                  - bg_vsno_cat[k](jnode, 0));
      if (dv > snowTouchedTol) {
        snowTouched[jnode * ncat_ + k] = 1;
      }
    }
  }

  applyThermoStage(pproc.fieldSet(), snowTouched);

  if (thermoParams.seedNewIce.value()) {
    const std::size_t fallbacks = seedNewIce(pproc.fieldSet(), newIce,
                                             donorCache, ice_lev, sno_lev);
    if (fallbacks > 0) {
      oops::Log::warning() << "PostProcessIce: noice->ice seed fell back to "
                           << "Tfrz at " << fallbacks << " cells (no donor "
                           << "with ice within seedSearchK)." << std::endl;
    }
  }
  // `pproc` (now carrying the mutated ice + thermo/pond fields) is written by
  // the caller via soca::State::write with a `cice template` key.
}

// -----------------------------------------------------------------------------

void PostProcessIce::print(std::ostream &) const {
}

// -----------------------------------------------------------------------------

double PostProcessIce::totalAice(const std::vector<atlas::array::ArrayView<double, 2>> & aiceCat,
                                 size_t jnode) const {
  double sum = 0.0;
  for (size_t icat = 0; icat < ncat_; ++icat) {
    sum += aiceCat[icat](jnode, 0);
  }
  return sum;
}

// -----------------------------------------------------------------------------

double PostProcessIce::meanHice(const std::vector<atlas::array::ArrayView<double, 2>> & viceCat,
                                double aice, size_t jnode) const {
  double vtot = 0.0;
  for (size_t icat = 0; icat < ncat_; ++icat) {
    vtot += viceCat[icat](jnode, 0);
  }
  return (aice > 0.0) ? std::max(0.0, vtot / aice) : 0.0;
}

// -----------------------------------------------------------------------------

double PostProcessIce::meanHsno(const std::vector<atlas::array::ArrayView<double, 2>> & vsnoCat,
                                double aice, size_t jnode) const {
  double vtot = 0.0;
  for (size_t icat = 0; icat < ncat_; ++icat) {
    vtot += vsnoCat[icat](jnode, 0);
  }
  return (aice > 0.0) ? std::max(0.0, vtot / aice) : 0.0;
}

// -----------------------------------------------------------------------------

void PostProcessIce::applyThermoStage(
    const atlas::FieldSet & fset,
    const std::vector<std::uint8_t> & snowTouched) const {
  const auto & thermoParams = params_.thermo.value();
  const bool updateSnowThermo = thermoParams.updateSnowThermo.value();
  const bool resetPonds       = thermoParams.resetPonds.value();
  if (!updateSnowThermo && !resetPonds) return;

  const std::size_t iceLev = params_.ice_lev.value();
  const std::size_t snoLev = params_.sno_lev.value();

  const double maxTsfc = thermoParams.maxTsfc.value();
  const double minTsfc = thermoParams.minTsfc.value();
  const double qsnoMax = icephysics::snowEnthalpy(maxTsfc);
  const double qsnoMin = icephysics::snowEnthalpy(minTsfc);
  const double qsnoLo  = std::min(qsnoMin, qsnoMax);
  const double qsnoHi  = std::max(qsnoMin, qsnoMax);

  auto aiceCat = ppiCatViews(fset, ncat_, "sea_ice_category", "_area_fraction");
  auto vsnoCat = ppiCatViews(fset, ncat_, "sea_ice_snow_category", "_volume");
  auto tsfcCat = ppiCatViews(fset, ncat_, "sea_ice_category", "_surface_temperature");
  auto apndCat = ppiCatViews(fset, ncat_, "sea_ice_category", "_pond_area_fraction");
  auto hpndCat = ppiCatViews(fset, ncat_, "sea_ice_category", "_pond_depth");
  auto ipndCat = ppiCatViews(fset, ncat_, "sea_ice_category", "_pond_lid_thickness");
  auto qiceCatLev = ppiCatLevViews(fset, ncat_, iceLev,
                                   "sea_ice_category", "_layer", "_enthalpy");
  auto siceCatLev = ppiCatLevViews(fset, ncat_, iceLev,
                                   "sea_ice_category", "_layer", "_ice_salinity");
  auto qsnoCatLev = ppiCatLevViews(fset, ncat_, snoLev,
                                   "sea_ice_snow_category", "_layer", "_enthalpy");

  const std::size_t field_size = aiceCat[0].shape(0);
  if (resetPonds) {
    ASSERT(snowTouched.size() == field_size * ncat_);
  }

  // Pre-compute layer-mean salinities so we can cap the surface ice layer.
  std::vector<double> sLayer(iceLev, 0.0);
  for (std::size_t l = 0; l < iceLev; ++l) {
    sLayer[l] = icephysics::siceLayerCice4(static_cast<int>(l)+1,
                                           static_cast<int>(iceLev));
  }
  const std::size_t lSurf = iceLev - 1;

  for (std::size_t jnode = 0; jnode < field_size; ++jnode) {
    for (std::size_t k = 0; k < ncat_; ++k) {
      const double aice = aiceCat[k](jnode, 0);
      if (aice <= 0.0) {
        // Cat lost its ice (ice2noice cell, rebin shuffle, or min-vice cleanup).
        // Clear stale thermo so cats without ice carry zeros rather than the
        // bg values they inherited from `pproc = restart`.
        tsfcCat[k](jnode, 0) = 0.0;
        for (std::size_t l = 0; l < snoLev; ++l) {
          qsnoCatLev[k][l](jnode, 0) = 0.0;
        }
        for (std::size_t l = 0; l < iceLev; ++l) {
          qiceCatLev[k][l](jnode, 0) = 0.0;
          siceCatLev[k][l](jnode, 0) = 0.0;
        }
        apndCat[k](jnode, 0) = 0.0;
        hpndCat[k](jnode, 0) = 0.0;
        ipndCat[k](jnode, 0) = 0.0;
        continue;
      }

      if (updateSnowThermo) {
        const double vsno = vsnoCat[k](jnode, 0);
        // Tsfcn from prior stages (preserved from background for ICE2ICE, or
        // donor-seeded by seedNewIce for NOICE2ICE) is the physically
        // meaningful quantity here. Rebuild qsno from Tsfcn rather than the
        // reverse, mirroring the Python reference (qsn_tsfc =
        // snowEnthalpy(tsf_new) in insert_iconc_ithkn_cice6_restart.py).
        // The opposite direction (back-deriving Tsfcn from a stale bg qsno
        // on cats that gained snow via redistribution) produced a +2.6 K
        // warm bias on the production grid.
        double & T = tsfcCat[k](jnode, 0);
        T = std::min(maxTsfc, std::max(minTsfc, T));
        if (vsno > 0.0) {
          const double qsfc = icephysics::snowEnthalpy(T);
          const double qClipped = std::min(qsnoHi, std::max(qsnoLo, qsfc));
          for (std::size_t l = 0; l < snoLev; ++l) {
            qsnoCatLev[k][l](jnode, 0) = qClipped;
          }
        } else {
          // No snow: zero qsno (matches Python; avoids stale-bg leakage).
          for (std::size_t l = 0; l < snoLev; ++l) {
            qsnoCatLev[k][l](jnode, 0) = 0.0;
          }
        }
        // Cap the surface ice layer enthalpy by iceEnthalpyBL99(Tsfcn, sice).
        if (T < 0.0) {
          const double sice = siceCatLev[k][lSurf](jnode, 0);
          const double sBL  = (sice > 0.0) ? sice : sLayer[lSurf];
          const double qCap = icephysics::iceEnthalpyBL99(T, sBL);
          double & qIce = qiceCatLev[k][lSurf](jnode, 0);
          qIce = std::min(qIce, qCap);
        }
      }

      if (resetPonds && snowTouched[jnode * ncat_ + k]) {
        // Match dd: only zero apnd/hpnd where snow was inserted. Never
        // touch ipnd (the refrozen-lid thickness is preserved from bg).
        apndCat[k](jnode, 0) = 0.0;
        hpndCat[k](jnode, 0) = 0.0;
      }
    }
  }
}

// -----------------------------------------------------------------------------

std::unordered_map<std::int64_t, PostProcessIce::CatRecord>
PostProcessIce::gatherDonorHalo(
    std::size_t K,
    const std::vector<atlas::array::ArrayView<double, 2>> & bg_aice_cat,
    const std::vector<atlas::array::ArrayView<double, 2>> & bg_vice_cat,
    const std::vector<atlas::array::ArrayView<double, 2>> & bg_vsno_cat,
    const std::vector<atlas::array::ArrayView<double, 2>> & bg_tsfc_cat,
    const std::vector<std::vector<atlas::array::ArrayView<double, 2>>> &
        bg_qice_cat_lev,
    const std::vector<std::vector<atlas::array::ArrayView<double, 2>>> &
        bg_sice_cat_lev,
    const std::vector<std::vector<atlas::array::ArrayView<double, 2>>> &
        bg_qsno_cat_lev,
    std::size_t ice_lev,
    std::size_t sno_lev) const {
  const auto & comm = geom_.getComm();
  const std::size_t nr = comm.size();
  const auto & fs   = geom_.functionSpace();
  const auto ghost  = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindex = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());

  // Phase A: per-cell KDTree query for K nearest neighbors. Accumulate the
  // union of returned gidx into a per-rank wanted set, dropping self.
  std::unordered_set<std::int64_t> wanted;
  for (std::size_t jnode = 0; jnode < ghost.size(); ++jnode) {
    if (ghost(jnode)) continue;
    if (gindex(jnode) <= 0) continue;
    const std::int64_t myG = static_cast<std::int64_t>(gindex(jnode));
    atlas::PointLonLat target(lonlat_(jnode, 0), lonlat_(jnode, 1));
    target.normalise();
    auto list = kdTree_.closestPoints(target, K);
    for (const auto & element : list) {
      const std::int64_t g = static_cast<std::int64_t>(element.payload());
      if (g != myG) wanted.insert(g);
    }
  }

  // Phase B: group wanted gidx by owner rank, allToAllv to send requests.
  std::vector<std::vector<atlas::gidx_t>> sendReq(nr);
  for (std::int64_t g : wanted) {
    auto it = gidxToOwnerRank_.find(g);
    if (it == gidxToOwnerRank_.end()) continue;  // shouldn't happen
    sendReq[it->second].push_back(static_cast<atlas::gidx_t>(g));
  }
  std::vector<std::vector<atlas::gidx_t>> recvReq(nr);
  comm.allToAll(sendReq, recvReq);

  // Phase C: for each requested gidx received, pack a flat reply record. Each
  // record packs in this fixed layout (recLen doubles per record):
  //   ncat aicen | ncat vicen | ncat vsnon | ncat Tsfcn |
  //   ncat*iceLev qice | ncat*iceLev sice | ncat*snoLev qsno | 1 mask
  const std::size_t recLen = ncat_ * 4
                           + ncat_ * ice_lev * 2
                           + ncat_ * sno_lev
                           + 1;
  std::vector<std::vector<double>> sendRep(nr);
  for (std::size_t r = 0; r < nr; ++r) {
    sendRep[r].resize(recvReq[r].size() * recLen, 0.0);
    for (std::size_t i = 0; i < recvReq[r].size(); ++i) {
      const std::int64_t g = static_cast<std::int64_t>(recvReq[r][i]);
      auto it = gidxToLocalJnode_.find(g);
      if (it == gidxToLocalJnode_.end()) {
        // Requester thinks we own this gidx but we don't — leave zeros.
        // Mask=0 (already zeroed) makes downstream skip it.
        continue;
      }
      const std::size_t jp = it->second;
      double * rec = sendRep[r].data() + i * recLen;
      // aicen, vicen, vsnon
      for (std::size_t k = 0; k < ncat_; ++k) {
        rec[k]              = bg_aice_cat[k](jp, 0);
        rec[ncat_ + k]      = bg_vice_cat[k](jp, 0);
        rec[2 * ncat_ + k]  = bg_vsno_cat[k](jp, 0);
      }
      // Tsfcn / per-layer qice/sice/qsno (per cat) from the background views.
      double * recTs = rec + 3 * ncat_;
      double * recQi = rec + 4 * ncat_;
      double * recSi = rec + 4 * ncat_ + ncat_ * ice_lev;
      double * recQs = rec + 4 * ncat_ + 2 * ncat_ * ice_lev;
      for (std::size_t k = 0; k < ncat_; ++k) {
        recTs[k] = bg_tsfc_cat[k](jp, 0);
        for (std::size_t l = 0; l < ice_lev; ++l) {
          recQi[k * ice_lev + l] = bg_qice_cat_lev[k][l](jp, 0);
          recSi[k * ice_lev + l] = bg_sice_cat_lev[k][l](jp, 0);
        }
        for (std::size_t l = 0; l < sno_lev; ++l) {
          recQs[k * sno_lev + l] = bg_qsno_cat_lev[k][l](jp, 0);
        }
      }
      // mask
      rec[recLen - 1] = mask_(jp, 0);
    }
  }
  std::vector<std::vector<double>> recvRep(nr);
  comm.allToAll(sendRep, recvRep);

  // Reassemble: each rank receives a flat list of records matching its
  // outgoing sendReq (per rank, in same order). Build the donor cache.
  std::unordered_map<std::int64_t, CatRecord> cache;
  cache.reserve(wanted.size());
  for (std::size_t r = 0; r < nr; ++r) {
    const std::size_t nrec = sendReq[r].size();
    ASSERT(recvRep[r].size() == nrec * recLen);
    for (std::size_t i = 0; i < nrec; ++i) {
      const std::int64_t g = static_cast<std::int64_t>(sendReq[r][i]);
      const double * rec = recvRep[r].data() + i * recLen;
      CatRecord cr;
      cr.aicen.assign(rec,                rec + ncat_);
      cr.vicen.assign(rec +     ncat_,    rec + 2 * ncat_);
      cr.vsnon.assign(rec + 2 * ncat_,    rec + 3 * ncat_);
      cr.Tsfcn.assign(rec + 3 * ncat_,    rec + 4 * ncat_);
      const double * pQi = rec + 4 * ncat_;
      const double * pSi = pQi + ncat_ * ice_lev;
      const double * pQs = pSi + ncat_ * ice_lev;
      cr.qice.assign(pQi, pQi + ncat_ * ice_lev);
      cr.sice.assign(pSi, pSi + ncat_ * ice_lev);
      cr.qsno.assign(pQs, pQs + ncat_ * sno_lev);
      cr.mask = rec[recLen - 1];
      cache.emplace(g, std::move(cr));
    }
  }
  return cache;
}

// -----------------------------------------------------------------------------

bool PostProcessIce::donorMeanIce(
    std::size_t jnode,
    std::size_t K,
    const std::unordered_map<std::int64_t, CatRecord> & donorCache,
    double & Tsfc_mean,
    double & hice_mean) const {
  const auto & fs   = geom_.functionSpace();
  const auto gindex = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());
  const std::int64_t myG = static_cast<std::int64_t>(gindex(jnode));

  atlas::PointLonLat target(lonlat_(jnode, 0), lonlat_(jnode, 1));
  target.normalise();
  auto list = kdTree_.closestPoints(target, K);

  for (const auto & element : list) {
    const std::int64_t g = static_cast<std::int64_t>(element.payload());
    if (g == myG) continue;
    auto it = donorCache.find(g);
    if (it == donorCache.end()) continue;
    const CatRecord & rec = it->second;
    if (rec.mask == 0.0) continue;
    double aiceTot   = 0.0;
    double tWeighted = 0.0;
    double vIceTot   = 0.0;
    for (std::size_t k = 0; k < ncat_; ++k) {
      if (rec.aicen[k] > 0.0) {
        aiceTot   += rec.aicen[k];
        tWeighted += rec.aicen[k] * rec.Tsfcn[k];
        vIceTot   += rec.vicen[k];
      }
    }
    if (aiceTot <= 0.0) continue;
    Tsfc_mean = tWeighted / aiceTot;
    hice_mean = vIceTot   / aiceTot;
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------

std::size_t PostProcessIce::seedNewIce(
    const atlas::FieldSet & fset,
    const NewIceMask & mask,
    const std::unordered_map<std::int64_t, CatRecord> & donorCache,
    std::size_t ice_lev,
    std::size_t sno_lev) const {
  const auto & thermoParams = params_.thermo.value();
  const std::size_t K = static_cast<std::size_t>(thermoParams.seedSearchK.value());
  const double maxTsfc = std::min(thermoParams.maxTsfc.value(), -1e-6);
  const double minTsfc = thermoParams.minTsfc.value();
  const double Tfrz    = icephysics::Constants::Tfrz;

  auto tsfcCat = ppiCatViews(fset, ncat_, "sea_ice_category", "_surface_temperature");
  auto qiceCatLev = ppiCatLevViews(fset, ncat_, ice_lev,
                                   "sea_ice_category", "_layer", "_enthalpy");
  auto siceCatLev = ppiCatLevViews(fset, ncat_, ice_lev,
                                   "sea_ice_category", "_layer", "_ice_salinity");
  auto qsnoCatLev = ppiCatLevViews(fset, ncat_, sno_lev,
                                   "sea_ice_snow_category", "_layer", "_enthalpy");

  const auto & fs   = geom_.functionSpace();
  const auto ghost  = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindex = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());

  std::size_t fallbackCount = 0;
  for (std::size_t jnode = 0; jnode < ghost.size(); ++jnode) {
    if (ghost(jnode)) continue;
    if (gindex(jnode) <= 0) continue;

    bool anyNew = false;
    for (std::size_t k = 0; k < ncat_; ++k) {
      if (mask.at(jnode, k)) { anyNew = true; break; }
    }
    if (!anyNew) continue;

    // Find the nearest-neighbor donor with any ice. Walk up to K candidates.
    atlas::PointLonLat target(lonlat_(jnode, 0), lonlat_(jnode, 1));
    target.normalise();
    auto list = kdTree_.closestPoints(target, K);
    const std::int64_t myG = static_cast<std::int64_t>(gindex(jnode));
    double donorTsfc = Tfrz;
    bool foundDonor = false;
    for (const auto & element : list) {
      const std::int64_t g = static_cast<std::int64_t>(element.payload());
      if (g == myG) continue;
      auto it = donorCache.find(g);
      if (it == donorCache.end()) continue;
      const CatRecord & rec = it->second;
      if (rec.mask == 0.0) continue;
      double aiceTot = 0.0;
      double tWeighted = 0.0;
      for (std::size_t k = 0; k < ncat_; ++k) {
        if (rec.aicen[k] > 0.0) {
          aiceTot   += rec.aicen[k];
          tWeighted += rec.aicen[k] * rec.Tsfcn[k];
        }
      }
      if (aiceTot <= 0.0) continue;
      donorTsfc = tWeighted / aiceTot;
      foundDonor = true;
      break;
    }
    if (!foundDonor) ++fallbackCount;

    const double Tseed = std::min(maxTsfc, std::max(minTsfc, donorTsfc));

    for (std::size_t k = 0; k < ncat_; ++k) {
      if (!mask.at(jnode, k)) continue;
      tsfcCat[k](jnode, 0) = Tseed;
      const double qsnSeed = icephysics::snowEnthalpy(Tseed);
      for (std::size_t l = 0; l < sno_lev; ++l) {
        qsnoCatLev[k][l](jnode, 0) = qsnSeed;
      }
      for (std::size_t l = 0; l < ice_lev; ++l) {
        const double sice_l = icephysics::siceLayerCice4(
            static_cast<int>(l) + 1, static_cast<int>(ice_lev));
        siceCatLev[k][l](jnode, 0) = sice_l;
        qiceCatLev[k][l](jnode, 0) =
            icephysics::iceEnthalpyBL99(Tseed, sice_l);
      }
    }
  }
  return fallbackCount;
}

// -----------------------------------------------------------------------------

}  // namespace soca
