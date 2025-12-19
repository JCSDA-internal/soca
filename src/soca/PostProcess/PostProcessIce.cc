/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/PostProcess/PostProcessIce.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/actions/BuildCellCentres.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
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
  // Build KDTree for local node indices
  const size_t npts = geom.functionSpace().size();
  std::vector<double> lons(npts);
  std::vector<double> lats(npts);
  std::vector<int> indices(npts);
  for (size_t jnode = 0; jnode < npts; ++jnode) {
    lons[jnode] = lonlat_(jnode, 0);
    lats[jnode] = lonlat_(jnode, 1);
    indices[jnode] = jnode;
  }
  kdTree_.build(lons, lats, indices);
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
  auto bg_tocn = atlas::array::make_view<double, 2>(bgfields.field("sea_water_potential_temperature"));
  auto bg_socn = atlas::array::make_view<double, 2>(bgfields.field("sea_water_salinity"));
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
    // Analysis has ice, but background has little ice: "shuffle"
      const size_t halo = std::max<size_t>(1, params_.shuffleStencilSize.value());
      const size_t stencilSize = (2*halo + 1) * (2*halo + 1);
      // search K nearest neighbors using KDTree
      atlas::PointLonLat target(lonlat_(jnode, 0), lonlat_(jnode, 1));
      target.normalise();
      auto list = kdTree_.closestPoints(target, stencilSize);
      // Choose donor among K nearest whose aice is closest to analysis
      double bestDiff = std::numeric_limits<double>::max();
      size_t bestJ = jnode;
      for (auto element : list) {
        const size_t jp = element.payload();
        if (mask_(jp, 0) == 0) continue; // skip land
        const double diff = std::abs(bg_aice(jp, 0) - a_aice(jnode, 0));
        if (diff < bestDiff) { bestDiff = diff; bestJ = jp; }
      }
      // Apply donor category distributions
      for (size_t icat = 0; icat < ncat_; ++icat) {
        new_aice_cat[icat](jnode, 0) = bg_aice_cat[icat](bestJ, 0);
        new_vice_cat[icat](jnode, 0) = bg_vice_cat[icat](bestJ, 0);
        new_vsno_cat[icat](jnode, 0) = bg_vsno_cat[icat](bestJ, 0);
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
  oops::Log::info() << " after pp restart: " << restart << std::endl;
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

}  // namespace soca
