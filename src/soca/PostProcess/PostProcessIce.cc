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
  auto bg_tocn = atlas::array::make_view<double, 2>(
      bgfields.field("sea_water_potential_temperature"));
  auto bg_socn = atlas::array::make_view<double, 2>(
      bgfields.field("sea_water_salinity"));
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

  // Analysis fields
  auto a_aice = atlas::array::make_view<double, 2>(anfields.field("sea_ice_area_fraction"));
  auto a_hice = atlas::array::make_view<double, 2>(anfields.field("sea_ice_thickness"));
  auto a_hsno = atlas::array::make_view<double, 2>(anfields.field("sea_ice_snow_thickness"));

  // parameters
  const auto & arctic = params_.arctic.value();
  const auto & antarctic = params_.antarctic.value();
  const double min_aice = params_.minAice.value();
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
  // CICE restart input: open + read thermo before the per-cell loop so the
  // shuffle path can copy donor thermo, and so the donor halo gather can pull
  // per-cat Tsfcn from each rank's local frame.
  const auto & ciceCfg = params_.ciceRestart.value();
  CiceRestartIO ciceIO(geom_, ciceCfg.input.value(), ciceCfg.output.value());
  auto frame = ciceIO.readThermo(ncat_, ice_lev, sno_lev);

  // ownedNodeOf[jnode] = owned-node index (frame slot) for `jnode`, or -1 for
  // ghost / non-owned cells. Mirrors CiceRestartIO::ownedGlobalIndices walk.
  std::vector<std::int64_t> ownedNodeOf(field_size, -1);
  {
    const auto ghost  = atlas::array::make_view<int, 1>(fs.ghost());
    const auto gindex = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());
    std::size_t on = 0;
    for (std::size_t jnode = 0; jnode < ghost.size(); ++jnode) {
      if (ghost(jnode)) continue;
      if (gindex(jnode) <= 0) continue;
      ownedNodeOf[jnode] = static_cast<std::int64_t>(on);
      ++on;
    }
    ASSERT(on == frame.nOwnedNodes);
  }

  // Sparse halo exchange for donor data (Phases A+B+C). K covers the maximum
  // of the shuffle stencil and the seed-new-ice search radius so both passes
  // can read from the same donorCache.
  const std::size_t halo     = std::max<std::size_t>(1,
                                 params_.shuffleStencilSize.value());
  const std::size_t stencilK = (2 * halo + 1) * (2 * halo + 1);
  const std::size_t seedK    = static_cast<std::size_t>(
      params_.thermo.value().seedSearchK.value());
  const std::size_t haloK    = std::max(stencilK, seedK);
  auto donorCache = gatherDonorHalo(haloK, bg_aice_cat, bg_vice_cat, bg_vsno_cat,
                                    frame, ownedNodeOf, ice_lev, sno_lev);

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
    const auto & reg = lonlat_(jnode, 1) > 0.0 ? arctic : antarctic;
    const bool adjust_sst = reg.adjustSST.value();
    const double max_dsst = reg.sstDiffMax.value();

    // Clamp bounds on analysis first.
    if (a_aice(jnode, 0) < min_aice) a_aice(jnode, 0) = 0.0;
    if (a_aice(jnode, 0) > 1.0)      a_aice(jnode, 0) = 1.0;
    if (a_hice(jnode, 0) < 0.0)      a_hice(jnode, 0) = 0.0;
    if (a_hsno(jnode, 0) < 0.0)      a_hsno(jnode, 0) = 0.0;
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
      double hice_seed = a_hice(jnode, 0);
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
      vtot_target = a_hice(jnode, 0) * ai_an;
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
      double hsn_target = a_hsno(jnode, 0);
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
  // Stage C: build NewIceMask, run applyThermoStage, optional seedNewIce,
  // write CICE output + flush thermo. With the shuffle path removed, every
  // (owned, k) cell that transitioned bg_aicen==0 -> new_aicen>0 needs thermo
  // seeding from a KDTree donor.
  NewIceMask newIce;
  newIce.ncat = ncat_;
  newIce.nOwnedNodes = frame.nOwnedNodes;
  newIce.data.assign(newIce.nOwnedNodes * ncat_, 0);
  for (std::size_t jnode = 0; jnode < field_size; ++jnode) {
    const std::int64_t on64 = ownedNodeOf[jnode];
    if (on64 < 0) continue;
    const std::size_t on = static_cast<std::size_t>(on64);
    for (std::size_t k = 0; k < ncat_; ++k) {
      const bool wasZero = (bg_aice_cat[k](jnode, 0) == 0.0);
      const bool nowPos  = (new_aice_cat[k](jnode, 0) > 0.0);
      if (wasZero && nowPos) {
        newIce.data[on * ncat_ + k] = 1;
      }
    }
  }

  const auto & thermoParams = params_.thermo.value();
  applyThermoStage(frame, pproc.fieldSet());

  if (thermoParams.seedNewIce.value()) {
    const std::size_t fallbacks = seedNewIce(frame, newIce, donorCache,
                                             ownedNodeOf, ice_lev, sno_lev);
    if (fallbacks > 0) {
      oops::Log::warning() << "PostProcessIce: noice->ice seed fell back to "
                           << "Tfrz at " << fallbacks << " cells (no donor "
                           << "with ice within seedSearchK)." << std::endl;
    }
  }

  ciceIO.write(pproc.fieldSet(), ncat_);
  ciceIO.flushThermo(frame);
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

void PostProcessIce::applyThermoStage(CiceRestartIO::ThermoFrame & frame,
                                      const atlas::FieldSet & fset) const {
  const auto & thermoParams = params_.thermo.value();
  const bool updateSnowThermo = thermoParams.updateSnowThermo.value();
  const bool resetPonds       = thermoParams.resetPonds.value();
  if (!updateSnowThermo && !resetPonds) return;

  ASSERT(frame.ncat == ncat_);
  const std::size_t iceLev = frame.iceLev;
  const std::size_t snoLev = frame.snoLev;

  const double maxTsfc = thermoParams.maxTsfc.value();
  const double minTsfc = thermoParams.minTsfc.value();
  const double qsnoMax = icephysics::snowEnthalpy(maxTsfc);
  const double qsnoMin = icephysics::snowEnthalpy(minTsfc);
  const double qsnoLo  = std::min(qsnoMin, qsnoMax);
  const double qsnoHi  = std::max(qsnoMin, qsnoMax);

  std::vector<atlas::array::ArrayView<double, 2>> aiceCat;
  std::vector<atlas::array::ArrayView<double, 2>> vsnoCat;
  aiceCat.reserve(ncat_);
  vsnoCat.reserve(ncat_);
  for (std::size_t icat = 0; icat < ncat_; ++icat) {
    const std::string aname = "sea_ice_category" + std::to_string(icat+1)
                              + "_area_fraction";
    const std::string sname = "sea_ice_snow_category" + std::to_string(icat+1)
                              + "_volume";
    aiceCat.push_back(atlas::array::make_view<double, 2>(fset.field(aname)));
    vsnoCat.push_back(atlas::array::make_view<double, 2>(fset.field(sname)));
  }

  // Pre-compute layer-mean salinities so we can cap the surface ice layer.
  std::vector<double> sLayer(iceLev, 0.0);
  for (std::size_t l = 0; l < iceLev; ++l) {
    sLayer[l] = icephysics::siceLayerCice4(static_cast<int>(l)+1,
                                           static_cast<int>(iceLev));
  }
  const std::size_t lSurf = iceLev - 1;

  const auto & fs    = geom_.functionSpace();
  const auto ghost   = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindex  = atlas::array::make_view<atlas::gidx_t, 1>(
                          fs.global_index());

  std::size_t ownedNode = 0;
  for (std::size_t jnode = 0; jnode < ghost.size(); ++jnode) {
    if (ghost(jnode)) continue;
    if (gindex(jnode) <= 0) continue;
    for (std::size_t k = 0; k < ncat_; ++k) {
      const double aice = aiceCat[k](jnode, 0);
      if (aice <= 0.0) continue;

      if (updateSnowThermo) {
        const double vsno = vsnoCat[k](jnode, 0);
        if (vsno > 0.0) {
          for (std::size_t l = 0; l < snoLev; ++l) {
            double & q = frame.at3(frame.qsno, snoLev, ownedNode, k, l);
            q = std::min(qsnoHi, std::max(qsnoLo, q));
          }
          // Use the surface snow layer to back out Tsfcn.
          const double qsfc = frame.at3(frame.qsno, snoLev, ownedNode, k, 0);
          const double tsfc = icephysics::snowEnthalpyToTsfc(qsfc);
          frame.at2(frame.Tsfcn, ownedNode, k) = std::min(maxTsfc, tsfc);
        }
        // Cap the surface ice layer enthalpy by iceEnthalpyBL99(Tsfcn, sice).
        const double tice = std::min(maxTsfc,
            frame.at2(frame.Tsfcn, ownedNode, k));
        if (tice < 0.0) {
          const double sice = frame.at3(frame.sice, iceLev,
                                         ownedNode, k, lSurf);
          const double sBL  = (sice > 0.0) ? sice : sLayer[lSurf];
          const double qCap = icephysics::iceEnthalpyBL99(tice, sBL);
          double & qIce = frame.at3(frame.qice, iceLev, ownedNode, k, lSurf);
          qIce = std::min(qIce, qCap);
        }
      }

      if (resetPonds) {
        frame.at2(frame.apnd, ownedNode, k) = 0.0;
        frame.at2(frame.hpnd, ownedNode, k) = 0.0;
        frame.at2(frame.ipnd, ownedNode, k) = 0.0;
      }
    }
    ++ownedNode;
  }
  ASSERT(ownedNode == frame.nOwnedNodes);
}

// -----------------------------------------------------------------------------

std::unordered_map<std::int64_t, PostProcessIce::CatRecord>
PostProcessIce::gatherDonorHalo(
    std::size_t K,
    const std::vector<atlas::array::ArrayView<double, 2>> & bg_aice_cat,
    const std::vector<atlas::array::ArrayView<double, 2>> & bg_vice_cat,
    const std::vector<atlas::array::ArrayView<double, 2>> & bg_vsno_cat,
    const CiceRestartIO::ThermoFrame & frame,
    const std::vector<std::int64_t> & ownedNodeOf,
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
      // Tsfcn (per cat) from local frame.
      const std::int64_t on = ownedNodeOf[jp];
      double * recTs = rec + 3 * ncat_;
      double * recQi = rec + 4 * ncat_;
      double * recSi = rec + 4 * ncat_ + ncat_ * ice_lev;
      double * recQs = rec + 4 * ncat_ + 2 * ncat_ * ice_lev;
      if (on >= 0) {
        const std::size_t onSt = static_cast<std::size_t>(on);
        for (std::size_t k = 0; k < ncat_; ++k) {
          recTs[k] = frame.at2(frame.Tsfcn, onSt, k);
          for (std::size_t l = 0; l < ice_lev; ++l) {
            recQi[k * ice_lev + l] = frame.at3(frame.qice, ice_lev, onSt, k, l);
            recSi[k * ice_lev + l] = frame.at3(frame.sice, ice_lev, onSt, k, l);
          }
          for (std::size_t l = 0; l < sno_lev; ++l) {
            recQs[k * sno_lev + l] = frame.at3(frame.qsno, sno_lev, onSt, k, l);
          }
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
    CiceRestartIO::ThermoFrame & frame,
    const NewIceMask & mask,
    const std::unordered_map<std::int64_t, CatRecord> & donorCache,
    const std::vector<std::int64_t> & ownedNodeOf,
    std::size_t ice_lev,
    std::size_t sno_lev) const {
  const auto & thermoParams = params_.thermo.value();
  const std::size_t K = static_cast<std::size_t>(thermoParams.seedSearchK.value());
  const double maxTsfc = std::min(thermoParams.maxTsfc.value(), -1e-6);
  const double minTsfc = thermoParams.minTsfc.value();
  const double Tfrz    = icephysics::Constants::Tfrz;

  const auto & fs   = geom_.functionSpace();
  const auto ghost  = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindex = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());

  std::size_t fallbackCount = 0;
  for (std::size_t jnode = 0; jnode < ghost.size(); ++jnode) {
    if (ghost(jnode)) continue;
    if (gindex(jnode) <= 0) continue;
    const std::int64_t on64 = ownedNodeOf[jnode];
    if (on64 < 0) continue;
    const std::size_t on = static_cast<std::size_t>(on64);

    bool anyNew = false;
    for (std::size_t k = 0; k < ncat_; ++k) {
      if (mask.at(on, k)) { anyNew = true; break; }
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
      if (!mask.at(on, k)) continue;
      frame.at2(frame.Tsfcn, on, k) = Tseed;
      const double qsnSeed = icephysics::snowEnthalpy(Tseed);
      for (std::size_t l = 0; l < sno_lev; ++l) {
        frame.at3(frame.qsno, sno_lev, on, k, l) = qsnSeed;
      }
      for (std::size_t l = 0; l < ice_lev; ++l) {
        const double sice_l = icephysics::siceLayerCice4(
            static_cast<int>(l) + 1, static_cast<int>(ice_lev));
        frame.at3(frame.sice, ice_lev, on, k, l) = sice_l;
        frame.at3(frame.qice, ice_lev, on, k, l) =
            icephysics::iceEnthalpyBL99(Tseed, sice_l);
      }
    }
  }
  return fallbackCount;
}

// -----------------------------------------------------------------------------

}  // namespace soca
