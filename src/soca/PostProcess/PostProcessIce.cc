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
  const bool snow_redistribute_area = snowParams.redistributeByArea.value();

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

  // Track which cells had their thermo seeded by the shuffle path so that the
  // post-loop seedNewIce can exclude them (donor consistency rule).
  std::vector<std::uint8_t> shuffleSeeded(field_size, 0);

  // Cache global_index / ghost views for fast lookup in the shuffle branch.
  const auto ghostView  = atlas::array::make_view<int, 1>(fs.ghost());
  const auto gindexView = atlas::array::make_view<atlas::gidx_t, 1>(fs.global_index());

  // Loop over all grid points to update ice fields
  for (size_t jnode = 0; jnode < field_size; ++jnode) {
    const auto & reg = lonlat_(jnode, 1) > 0.0 ? arctic : antarctic;
    const double edge = reg.edge.value();
    const bool do_shuffle = reg.shuffle.value();
    const bool do_rescale = reg.rescale.value().rescale.value();
    const double min_hice = reg.rescale.value().min_hice.value();
    const double min_hsno = reg.rescale.value().min_hsno.value();
    const bool adjust_sst = reg.adjustSST.value();
    const double max_dsst = reg.sstDiffMax.value();

    // Clamp bounds on analysis first
    // ice concentration bounds
    if (a_aice(jnode, 0) < min_aice) a_aice(jnode, 0) = 0.0;
    if (a_aice(jnode, 0) > 1.0)      a_aice(jnode, 0) = 1.0;
    // non-negative thicknesses
    if (a_hice(jnode, 0) < 0.0) a_hice(jnode, 0) = 0.0;
    if (a_hsno(jnode, 0) < 0.0) a_hsno(jnode, 0) = 0.0;

    // Zero out ice if:
    // - analysis has no ice
    // - it's a land point
    if ((mask_(jnode, 0) == 0.0) || (a_aice(jnode, 0) <= 0.0)) {
      for (size_t icat = 0; icat < ncat_; ++icat) {
        new_aice_cat[icat](jnode, 0) = 0.0;
        new_vice_cat[icat](jnode, 0) = 0.0;
        new_vsno_cat[icat](jnode, 0) = 0.0;
      }
      // TODO(AS): adjust SST if needed
    } else if ((bg_aice(jnode, 0) <= edge) && do_shuffle) {
      // Analysis has ice, but background has little ice: "shuffle".
      // Pick donor among the K nearest neighbors whose total aice is closest
      // to the analysis aice. Donor data comes from donorCache (sparse halo
      // exchange) so cross-rank donors work correctly.
      const std::size_t shuffleK = (2 * halo + 1) * (2 * halo + 1);
      atlas::PointLonLat target(lonlat_(jnode, 0), lonlat_(jnode, 1));
      target.normalise();
      auto list = kdTree_.closestPoints(target, shuffleK);
      const std::int64_t myG = (gindexView(jnode) > 0)
          ? static_cast<std::int64_t>(gindexView(jnode)) : -1;
      double bestDiff = std::numeric_limits<double>::max();
      const CatRecord * bestRec = nullptr;
      for (const auto & element : list) {
        const std::int64_t g = static_cast<std::int64_t>(element.payload());
        if (g == myG) continue;
        auto it = donorCache.find(g);
        if (it == donorCache.end()) continue;
        const CatRecord & rec = it->second;
        if (rec.mask == 0.0) continue;  // skip land
        const double donorAice = std::accumulate(rec.aicen.begin(),
                                                 rec.aicen.end(), 0.0);
        const double diff = std::abs(donorAice - a_aice(jnode, 0));
        if (diff < bestDiff) { bestDiff = diff; bestRec = &rec; }
      }
      if (bestRec) {
        // Copy donor ice fields.
        for (std::size_t k = 0; k < ncat_; ++k) {
          new_aice_cat[k](jnode, 0) = bestRec->aicen[k];
          new_vice_cat[k](jnode, 0) = bestRec->vicen[k];
          new_vsno_cat[k](jnode, 0) = bestRec->vsnon[k];
        }
        // Copy donor thermo into the local ThermoFrame slot. Donor consistency
        // rule: this cell's thermo is now authoritative; seedNewIce will skip
        // it (shuffleSeeded[jnode] = true).
        const std::int64_t on64 = ownedNodeOf[jnode];
        if (on64 >= 0) {
          const std::size_t on = static_cast<std::size_t>(on64);
          for (std::size_t k = 0; k < ncat_; ++k) {
            frame.at2(frame.Tsfcn, on, k) = bestRec->Tsfcn[k];
            for (std::size_t l = 0; l < ice_lev; ++l) {
              frame.at3(frame.qice, ice_lev, on, k, l) =
                  bestRec->qice[k * ice_lev + l];
              frame.at3(frame.sice, ice_lev, on, k, l) =
                  bestRec->sice[k * ice_lev + l];
            }
            for (std::size_t l = 0; l < sno_lev; ++l) {
              frame.at3(frame.qsno, sno_lev, on, k, l) =
                  bestRec->qsno[k * sno_lev + l];
            }
          }
          shuffleSeeded[jnode] = 1;
        }
      }
    }
    // Rescale ice distribution if needed
    if ((mask_(jnode, 0) > 0.0) && (a_aice(jnode, 0) > 0.0) && do_rescale) {
      // Rescale ice concentration
      double alpha = a_aice(jnode, 0) / bg_aice(jnode, 0);
      new_aice(jnode, 0) = a_aice(jnode, 0);
      for (size_t icat = 0; icat < ncat_; ++icat) {
        new_aice_cat[icat](jnode, 0) = alpha * bg_aice_cat[icat](jnode, 0);
        new_vice_cat[icat](jnode, 0) = alpha * bg_vice_cat[icat](jnode, 0);
        new_vsno_cat[icat](jnode, 0) = alpha * bg_vsno_cat[icat](jnode, 0);
      }
      // Compute background mean thicknesses
      double hice_bg = meanHice(new_vice_cat, new_aice(jnode, 0), jnode);
      double hsno_bg = meanHsno(new_vsno_cat, new_aice(jnode, 0), jnode);
      // Rescale thicknesses to match analysis
      if (hice_bg > min_hice) {
        alpha = a_hice(jnode, 0) / hice_bg;
        // Rescale category volumes
        for (size_t icat = 0; icat < ncat_; ++icat) {
          new_vice_cat[icat](jnode, 0) = alpha * new_vice_cat[icat](jnode, 0);
        }
      }
      if (hsno_bg > min_hsno) {
        alpha = a_hsno(jnode, 0) / hsno_bg;
        // Rescale category snow volumes
        for (size_t icat = 0; icat < ncat_; ++icat) {
          new_vsno_cat[icat](jnode, 0) = alpha * new_vsno_cat[icat](jnode, 0);
        }
      }
    }
    // Post-rescale snow refinement: optional area-weighted redistribution and
    // a minimum-thickness floor. Both are no-ops by default to preserve parity
    // with the legacy Fortran path.
    if ((mask_(jnode, 0) > 0.0) && (a_aice(jnode, 0) > 0.0)
        && (snow_redistribute_area || hsnow_min > 0.0)) {
      // Snow volume currently in this cell; aice already matches analysis.
      double aice_now = 0.0;
      double vsno_sum = 0.0;
      for (size_t icat = 0; icat < ncat_; ++icat) {
        aice_now += new_aice_cat[icat](jnode, 0);
        vsno_sum += new_vsno_cat[icat](jnode, 0);
      }
      if (aice_now > 0.0) {
        if (snow_redistribute_area) {
          // Distribute total snow volume by aicen weight: vsnon[k] = vsno_sum
          // * aicen[k] / aice_now. Conserves Σ vsnon and keeps mean hsno =
          // vsno_sum / aice_now.
          for (size_t icat = 0; icat < ncat_; ++icat) {
            new_vsno_cat[icat](jnode, 0) =
                vsno_sum * (new_aice_cat[icat](jnode, 0) / aice_now);
          }
        }
        if (hsnow_min > 0.0) {
          // Drop snow on any cat whose mean snow thickness is below the floor.
          for (size_t icat = 0; icat < ncat_; ++icat) {
            const double a = new_aice_cat[icat](jnode, 0);
            if (a > 0.0 && (new_vsno_cat[icat](jnode, 0) / a) < hsnow_min) {
              new_vsno_cat[icat](jnode, 0) = 0.0;
            }
          }
        }
      }
    }
    // Re-bin the ice thickness distribution so each category's mean thickness
    // is inside its bin. Keeps the per-cell totals (sum aicen, sum vicen) the
    // same; only redistributes mass across categories to satisfy the bounds.
    if (do_rebin && (mask_(jnode, 0) > 0.0)) {
      double aice_sum = 0.0;
      double vice_sum = 0.0;
      for (size_t icat = 0; icat < ncat_; ++icat) {
        rebin_aicen[icat] = new_aice_cat[icat](jnode, 0);
        rebin_vicen[icat] = new_vice_cat[icat](jnode, 0);
        aice_sum += rebin_aicen[icat];
        vice_sum += rebin_vicen[icat];
      }
      if (aice_sum > 0.0 && vice_sum > 0.0) {
        const bool ok = icephysics::adjustThicknessCategories(
            rebin_aicen, rebin_vicen, vice_sum, hicat, dhi_min);
        if (ok) {
          for (size_t icat = 0; icat < ncat_; ++icat) {
            new_vice_cat[icat](jnode, 0) = rebin_vicen[icat];
          }
        } else {
          ++rebin_failures;
        }
      }
    }

    // Enforce hydrostatic snow/ice balance: rho_ice*hi + rho_snow*hs <=
    // rho_ocean*(hi+hs) per category. Snow is first redistributed across
    // cats; if any cat is still flooded, ice volume grows to lift the snow-
    // ice interface back to sea level. Per-cell, no neighbour search needed.
    if (do_freeboard && (mask_(jnode, 0) > 0.0)) {
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

    // Zero-out ice categories with tiny volumes
    for (size_t icat = 0; icat < ncat_; ++icat) {
      if (new_aice_cat[icat](jnode, 0) > 0.0 && new_vice_cat[icat](jnode, 0) < min_vice) {
        new_aice_cat[icat](jnode, 0) = 0.0;
        new_vice_cat[icat](jnode, 0) = 0.0;
        new_vsno_cat[icat](jnode, 0) = 0.0;
      }
    }
    // Recompute total ice concentration and mean thicknesses
    new_aice(jnode, 0) = totalAice(new_aice_cat, jnode);
    new_hice(jnode, 0) = meanHice(new_vice_cat, new_aice(jnode, 0), jnode);
    new_hsno(jnode, 0) = meanHsno(new_vsno_cat, new_aice(jnode, 0), jnode);
  }
  if (do_rebin && rebin_failures > 0) {
    oops::Log::warning() << "PostProcessIce: ITD rebin failed at "
                         << rebin_failures << " cells (target outside the "
                         << "feasible envelope); left untouched." << std::endl;
  }
  if (do_freeboard && freeboard_failures > 0) {
    oops::Log::warning() << "PostProcessIce: freeboard enforcement failed at "
                         << freeboard_failures << " cells; left untouched."
                         << std::endl;
  }
  oops::Log::info() << " after pp restart: " << restart << std::endl;

  // ---------------------------------------------------------------------------
  // Stage C: build NewIceMask, run applyThermoStage, optional seedNewIce,
  // write CICE output + flush thermo. Donor consistency rule: cells whose
  // thermo was already seeded by the shuffle path are excluded from the
  // newIce mask so seedNewIce doesn't re-seed them with a different donor.
  NewIceMask newIce;
  newIce.ncat = ncat_;
  newIce.nOwnedNodes = frame.nOwnedNodes;
  newIce.data.assign(newIce.nOwnedNodes * ncat_, 0);
  for (std::size_t jnode = 0; jnode < field_size; ++jnode) {
    const std::int64_t on64 = ownedNodeOf[jnode];
    if (on64 < 0) continue;
    if (shuffleSeeded[jnode]) continue;
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
