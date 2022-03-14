/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Interpolator/LocalUnstructuredInterpolator.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

#include "oops/util/abor1_cpp.h"

namespace soca {

// ------------------------------------------------------------------------------

const std::vector<char> grids{ 'h', 'u', 'v'};

LocalUnstructuredInterpolator::
  LocalUnstructuredInterpolator(const eckit::Configuration & config, const Geometry & geom,
                                const std::vector<double> & latlon_out)
  : geom_(new Geometry(geom)) {

  // fspace_(grid.atlasFunctionSpaceIncludingHalo())

  oops::Log::trace() << "LocalUnstructuredInterpolator::LocalUnstructuredInterpolator start" << std::endl;

  // create interpolator for each mask/staggering combination
  // TODO delay until needed to prevent u/v from being created if not used
  int idx = 0;
  for (char g: grids ) {
    // unmasked
    std::vector<double> lats_in;
    std::vector<double> lons_in;
    geom.latlon(lats_in, lons_in, true, g, false);
    interp_[idx++] = std::make_unique<oops::InterpolatorUnstructured>(config, lats_in, lons_in, latlon_out);

    // masked
    lats_in.clear();
    lons_in.clear();
    geom.latlon(lats_in, lons_in, true, g, true);
    interp_[idx++] = std::make_unique<oops::InterpolatorUnstructured>(config, lats_in, lons_in, latlon_out);
  }

  oops::Log::trace() << "LocalUnstructuredInterpolator::LocalUnstructuredInterpolator done" << std::endl;
  }

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::apply(
    const oops::Variables & vars, const State & xx, std::vector<double> & locvals )
    const {

  oops::Log::trace() << "LocalUnstructuredInterpolator::apply STATE start" << std::endl;
  locvals.clear();

  for (int i =0; i < vars.size(); i++) {
    // determine which interpolator to use
    char grid;
    bool masked;
    geom_->getVarGrid(vars[i], grid, masked);
    int interp_idx = -1;
    for (int j=0; j < grids.size(); j++) {
      if (grids[j] == grid) interp_idx = j*2;
    }
    if (masked) interp_idx += 1;

    // get a single variable
    oops::Variables var;
    var.push_back(vars[i]);
    atlas::FieldSet fset;
    xx.getFieldSet(var, fset);

    // interpolate
    std::vector<double> var_locvals;
    interp_[interp_idx]->apply(fset, var_locvals);
    locvals.insert(locvals.end(), var_locvals.begin(), var_locvals.end());
  }

  oops::Log::trace() << "LocalUnstructuredInterpolator::apply STATE done" << std::endl;

  }

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::
  apply(const oops::Variables & vars, const Increment & dx, std::vector<double> & locvals) const {
    util::abor1_cpp("not implemented: LocalUnstructuredInterpolator::apply(2)");
  }

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::
  applyAD(const oops::Variables & vars, Increment & dx, const std::vector<double> & locvals) const {
    util::abor1_cpp("not implemented: LocalUnstructuredInterpolator::applyAD()");
  }

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::print(std::ostream & os) const {

}

}