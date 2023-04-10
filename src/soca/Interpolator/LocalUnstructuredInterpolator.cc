/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <numeric>

#include "atlas/field.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/abor1_cpp.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Interpolator/UnstructuredInterpolator.h"
#include "soca/Interpolator/LocalUnstructuredInterpolator.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"



namespace soca {

// ------------------------------------------------------------------------------

const std::vector<char> grids{ 'h', 'u', 'v'};

LocalUnstructuredInterpolator::
  LocalUnstructuredInterpolator(const eckit::Configuration & config, const Geometry & geom,
                                const std::vector<double> & lats_out,
                                const std::vector<double> & lons_out)
  : geom_(geom) {
  // create interpolators (TODO, postpone to usage, or configure which interpolators
  // are created, since we likely don't need them all?)
  for (auto grid : grids) {
    util::Timer timer("soca::LocalUnstructuredInterpolator", "getInterpolator.build");
    int interp_idx = -1;
    for (int j=0; j < grids.size(); j++) {
      if (grids[j] == grid) interp_idx = j*2;
    }

    bool masked = false;
    interp_[interp_idx] = std::make_shared<UnstructuredInterpolator>(
        config, geom_, grid, masked, lats_out, lons_out);

    masked = true;
    interp_idx += 1;
    interp_[interp_idx] = std::make_shared<UnstructuredInterpolator>(
        config, geom_, grid, masked, lats_out, lons_out);
  }
}

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::
apply(const oops::Variables & vars, const atlas::FieldSet & fset, const std::vector<bool> & mask,
       std::vector<double> & locvals) const {
  auto vals = locvals.begin();
  for (int i =0; i < vars.size(); i++) {
    auto interpolator = getInterpolator(vars[i]);

    // get a single variable
    oops::Variables var;
    var.push_back(vars[i]);

    // interpolate
    interpolator->apply(var, fset, mask, vals);
  }
}

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::
applyAD(const oops::Variables & vars, atlas::FieldSet & fset, const std::vector<bool> & mask,
        const std::vector<double> & locvals) const {
  auto vals = locvals.begin();
  for (int i =0; i < vars.size(); i++) {
    auto interpolator = getInterpolator(vars[i]);

    // get a single variable
    oops::Variables var;
    var.push_back(vars[i]);

    // interpolate
    interpolator->applyAD(var, fset, mask, vals);
  }
}

// ------------------------------------------------------------------------------
const std::shared_ptr<UnstructuredInterpolator>
LocalUnstructuredInterpolator::getInterpolator(const std::string &var) const {
  // determine which interpolator to use
  char grid;
  bool masked;
  geom_.getVarGrid(var, grid, masked);
  int interp_idx = -1;
  for (int j=0; j < grids.size(); j++) {
    if (grids[j] == grid) interp_idx = j*2;
  }
  if (masked) interp_idx += 1;
  ASSERT(interp_idx < 6);

  return interp_[interp_idx];
}

void LocalUnstructuredInterpolator::print(std::ostream & os) const {
}

}  // namespace soca
