/*
* (C) Copyright 2026 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"

#include "oops/util/Logger.h"

#include "oops/base/Variables.h"

#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

#include "soca/Utils/incrqc/include/soca_diagb_utils.h"
#include "soca/Utils/incrqc/include/soca_physlight_utils.h"
#include "soca/Utils/incrqc/include/soca_incr_qc_utils.h"

namespace soca {
namespace incrqc {

/**
 * @brief Throws a helpful error if a field required by a QC stage is missing.
 *
 * @param fs The field set that should contain the field.
 * @param name The name of the required field.
 * @param fsName Human readable name of the field set ("increment", "background", ...).
 * @param stage The QC stage that requires the field, named as its configuration key.
 */
inline void requireField(const atlas::FieldSet& fs,
                         const std::string& name,
                         const std::string& fsName,
                         const std::string& stage) {
  if (!fs.has(name)) {
    oops::Log::error() << "Increment QC: the '" << stage << "' stage requires '" << name
                       << "' in the " << fsName << ", which does not contain it." << std::endl;
    throw eckit::UserError("Increment QC: '" + stage + "' requires '" + name +
                           "' in the " + fsName, Here());
  }
}

/**
 * @brief Throws a helpful error if a variable required by a QC stage is missing.
 *
 * @param vars The variables of the state or increment a QC stage will be applied to.
 * @param names The names of the required variables.
 * @param varsName Human readable name of what the variables belong to ("increment", ...).
 * @param stage The QC stage that requires the variables, named as its configuration key.
 */
inline void requireVariables(const oops::Variables& vars,
                             const std::vector<std::string>& names,
                             const std::string& varsName,
                             const std::string& stage) {
  for (const std::string & name : names) {
    if (!vars.has(name)) {
      oops::Log::error() << "Increment QC: the '" << stage << "' stage requires '" << name
                         << "' in the " << varsName << " variables, which are " << vars
                         << std::endl;
      throw eckit::UserError("Increment QC: '" + stage + "' requires '" + name +
                             "' in the " + varsName + " variables", Here());
    }
  }
}

/**
 * @brief Reads the per-variable analysis bounds from the "state bounds" configuration.
 *
 * Any variable can be given bounds, as a [minimum, maximum] pair. Variables that are
 * bounded but absent from the increment are reported and ignored, since a variable name
 * that does not match anything is much more likely to be a typo than to be intentional.
 *
 * @param config The QC configuration, containing a "state bounds" sub-configuration.
 * @param dxFs The increment field set, used to report bounds that cannot be applied.
 * @return Map of variable name to its {minimum, maximum} analysis bounds.
 */
inline std::unordered_map<std::string, std::pair<double, double>> getStateBounds(
    const eckit::Configuration& config,
    const atlas::FieldSet& dxFs) {
  const eckit::LocalConfiguration boundsConfig(config, "state bounds");
  std::unordered_map<std::string, std::pair<double, double>> stateBounds;

  for (const std::string & var : boundsConfig.keys()) {
    const std::vector<double> bounds = boundsConfig.getDoubleVector(var);
    if (bounds.size() != 2 || bounds[0] > bounds[1]) {
      oops::Log::error() << "Increment QC: 'state bounds." << var << "' must be a "
                         << "[minimum, maximum] pair with minimum <= maximum." << std::endl;
      throw eckit::UserError("Increment QC: invalid state bounds for '" + var + "'", Here());
    }
    if (!dxFs.has(var)) {
      oops::Log::warning() << "Increment QC: bounds given for '" << var << "', which is not "
                           << "an increment variable. The bounds are ignored." << std::endl;
      continue;
    }
    stateBounds[var] = {bounds[0], bounds[1]};
    oops::Log::info() << "QC: bounding the analysis of " << var
                      << " to [" << bounds[0] << ", " << bounds[1] << "]" << std::endl;
  }
  return stateBounds;
}

/**
 * @brief Reduces the per-variable clamp counts over all MPI tasks and logs them.
 *
 * The counts make it visible whether the bounds of each variable actually bit, which is
 * otherwise indistinguishable from the variable not being bounded at all.
 *
 * @param clampCounts Task-local {variable, {number clamped to minimum, to maximum}} counts.
 * @param comm The geometry communicator the counts are summed over.
 */
inline void logClampCounts(
    const std::map<std::string, std::pair<size_t, size_t>>& clampCounts,
    const eckit::mpi::Comm& comm) {
  for (const auto & entry : clampCounts) {
    std::vector<size_t> counts = {entry.second.first, entry.second.second};
    comm.allReduceInPlace(counts.begin(), counts.end(), eckit::mpi::sum());
    oops::Log::info() << "QC: state bounds on " << entry.first << ": clamped "
                      << counts[0] << " points to the minimum and "
                      << counts[1] << " points to the maximum" << std::endl;
  }
}

/**
 * @brief Quality control for increments: ensures that the analysis (xb + dx) remains within physical bounds.
 *
 * Every QC stage is optional and runs only if its own configuration key is present:
 * - "steric variable change": replace the ssh increment with a steric height increment, then
 *   limit it to "absolute steric increment max" by rescaling the T/S increments.
 * - "coastal increment filter": taper the T/S increments to zero near coastlines.
 * - "stability check": damp increments that destabilize the water column.
 * - "state bounds": keep the analysis of each listed variable within its bounds.
 *
 * @param xb The background state.
 * @param dx The increment to QC. Will be modified in place.
 * @param config The configuration containing bounds information.
 */
inline void qcIncrement(const soca::State& xb,
                 soca::Increment& dx,
                 const eckit::Configuration& config,
                 const soca::Geometry& geom) {
  oops::Log::info() << "==========================================" << std::endl;
  oops::Log::info() << "======      Quality control on increment" << std::endl;

  // Each stage is gated on the presence of its own configuration key
  const bool doSteric    = config.has("steric variable change");
  const bool doCoastal   = config.has("coastal increment filter");
  const bool doStability = config.has("stability check");
  const bool doBounds    = config.has("state bounds");

  oops::Log::info() << "QC: steric height adjustment: " << (doSteric ? "on" : "off") << std::endl;
  oops::Log::info() << "QC: coastal increment filter: " << (doCoastal ? "on" : "off") << std::endl;
  oops::Log::info() << "QC: water column stability check: "
                    << (doStability ? "on" : "off") << std::endl;
  oops::Log::info() << "QC: state bounds: " << (doBounds ? "on" : "off") << std::endl;

  if (!doSteric && !doCoastal && !doStability && !doBounds) {
    oops::Log::warning() << "Increment QC: no QC stage is configured, the increment is left "
                         << "untouched. Expected at least one of 'steric variable change', "
                         << "'coastal increment filter', 'stability check' or 'state bounds'."
                         << std::endl;
    return;
  }

  // Check that each enabled stage has the variables it needs, before anything is modified.
  // Bathymetry, and hence the layer thickness it is derived from, is needed by every stage
  // for the land mask.
  const std::vector<std::string> tempSalt = {"sea_water_potential_temperature",
                                             "sea_water_salinity"};
  requireVariables(xb.variables(), {"sea_water_cell_thickness"}, "background", "increment QC");
  if (doSteric) {
    requireVariables(dx.variables(), {"sea_water_potential_temperature", "sea_water_salinity",
                                      "sea_surface_height_above_geoid"},
                     "increment", "steric variable change");
    requireVariables(xb.variables(), tempSalt, "background", "steric variable change");
  }
  if (doCoastal) {
    requireVariables(dx.variables(), tempSalt, "increment", "coastal increment filter");
  }
  if (doStability) {
    requireVariables(dx.variables(), tempSalt, "increment", "stability check");
    requireVariables(xb.variables(), tempSalt, "background", "stability check");
  }

  // Replace the ssh increment with a steric height increment. This has to happen before
  // dx is copied into a field set below, since it modifies dx through a variable change.
  if (doSteric) {
    eckit::LocalConfiguration lvcConfig(config, "steric variable change");
    soca::utils::computeStericHeightIncrement(geom, dx, lvcConfig, xb, dx.variables());
  }

  atlas::FieldSet xbFs, dxFs;
  xb.toFieldSet(xbFs);
  dx.toFieldSet(dxFs);

  // Compute ocean depth and bathymetry
  auto viewHocn = atlas::array::make_view<double, 2>(xbFs["sea_water_cell_thickness"]);
  atlas::array::ArrayT<double> depth(viewHocn.shape(0), viewHocn.shape(1));
  auto viewDepth = atlas::array::make_view<double, 2>(depth);
  atlas::array::ArrayT<double> bathy(viewHocn.shape(0), 1);
  auto viewBathy = atlas::array::make_view<double, 2>(bathy);
  soca::utils::computeDepthAndBathymetry(viewHocn, viewDepth, viewBathy);

  // Get ghost nodes and lon/lat coordinates
  soca::diagb::utils::MeshBundle meshConn = soca::diagb::utils::buildMeshConnectivity(geom);
  const auto ghostView = meshConn.ghostView;
  const auto & lonlat = meshConn.lonlat;

  // Update halos for increment fields. Fields absent from the increment are skipped:
  // the stages that need them are off, this was checked above.
  const std::vector<std::string> fieldsToExchange = {
    "sea_water_potential_temperature",
    "sea_water_salinity",
    "sea_surface_height_above_geoid",
  };
  for (const auto& field : fieldsToExchange) {
    if (dxFs.has(field)) meshConn.nodeColumns.haloExchange(dxFs[field]);
  }

  // Update halos for background fields
  if (xbFs.has("sea_water_potential_temperature")) {
    meshConn.nodeColumns.haloExchange(xbFs["sea_water_potential_temperature"]);
  }
  if (xbFs.has("sea_water_salinity")) {
    meshConn.nodeColumns.haloExchange(xbFs["sea_water_salinity"]);
  }
  meshConn.nodeColumns.haloExchange(xbFs["sea_water_cell_thickness"]);

  // Coastal increment filter: taper T/S increments to zero near coastlines
  if (doCoastal) {
    const double distMin = config.getDouble("coastal increment filter.min distance");
    const double distMax = config.getDouble("coastal increment filter.max distance");
    // distance_from_coast lives in the geometry fieldset, not the state
    requireField(geom.fields(), "distance_from_coast", "geometry",
                 "coastal increment filter");
    auto viewDist = atlas::array::make_view<const double, 2>(
        geom.fields()["distance_from_coast"]);
    applyCoastalIncrementFilter(dxFs, viewDist, ghostView, viewBathy,
                                distMin, distMax);
  }

  // Stability checks
  if (doStability) {
    int niterations = config.getInt("stability check.increment stability iterations", 10);
    int nSmoothingIterations = config.getInt("stability check.increment smoothing iterations", 30);
    const double rhoMinGrad = config.getDouble("stability check.min stable density gradient", 1e-4);
    auto viewTempBkg = atlas::array::make_view<double, 2>(xbFs["sea_water_potential_temperature"]);
    auto viewSaltBkg = atlas::array::make_view<double, 2>(xbFs["sea_water_salinity"]);
    applyWaterColumnStabilityCheck(dxFs,
                                   viewTempBkg, viewSaltBkg,
                                   viewHocn, viewDepth, lonlat,
                                   niterations, rhoMinGrad, nSmoothingIterations,
                                   viewBathy, meshConn);
  }

  // Apply steric height constraint
  // Limit the SSH increment to deltaSshMax by rescaling T/S increments if necessary
  // Replaces the ssh increment with the final steric height increment
  if (doSteric) {
    const double deltaSshMax = config.getDouble("absolute steric increment max", 10.0);
    oops::Log::info() << "QC: max steric height increment: " << deltaSshMax << std::endl;

    auto viewTempIncr =
        atlas::array::make_view<double, 2>(dxFs["sea_water_potential_temperature"]);
    auto viewSaltIncr = atlas::array::make_view<double, 2>(dxFs["sea_water_salinity"]);
    auto viewSshIncr =
        atlas::array::make_view<double, 2>(dxFs["sea_surface_height_above_geoid"]);
    auto viewTempBkg = atlas::array::make_view<double, 2>(xbFs["sea_water_potential_temperature"]);
    auto viewSaltBkg = atlas::array::make_view<double, 2>(xbFs["sea_water_salinity"]);

    for (atlas::idx_t jnode = 0; jnode < viewTempIncr.shape(0); ++jnode) {
      // Skip ghost and land nodes
      if (ghostView(jnode) > 0) continue;
      if (viewBathy(jnode, 0) <= 0.0) continue;

      // Limit the steric height incrememnt to deltaSshMax
      applyStericHeightConstraint(jnode, viewTempIncr, viewSaltIncr, viewSshIncr,
                              viewTempBkg, viewSaltBkg, viewHocn, deltaSshMax);
    }
  }

  // Brute force bounds check
  if (doBounds) {
    const std::unordered_map<std::string, std::pair<double, double>> stateBounds =
        getStateBounds(config, dxFs);
    std::map<std::string, std::pair<size_t, size_t>> clampCounts =
        applyBruteForceBoundsCheck(dxFs, xbFs, ghostView, viewBathy, stateBounds);
    logClampCounts(clampCounts, geom.getComm());
  }

  dx.fromFieldSet(dxFs);
  oops::Log::info() << "======      Finished quality control on increment" << std::endl;
  }

}  // namespace incrqc
}  // namespace soca
