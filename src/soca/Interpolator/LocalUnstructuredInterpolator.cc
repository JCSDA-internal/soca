/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <numeric>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Interpolator/UnstructuredInterpolator.h"
#include "soca/Interpolator/LocalUnstructuredInterpolator.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

#include "oops/util/abor1_cpp.h"

namespace soca {

// ------------------------------------------------------------------------------

const std::vector<char> grids{ 'h', 'u', 'v'};

LocalUnstructuredInterpolator::
  LocalUnstructuredInterpolator(const eckit::Configuration & config, const Geometry & geom,
                                const std::vector<double> & lats_out,
                                const std::vector<double> & lons_out)
  : config_(config), lats_out_(lats_out), lons_out_(lons_out), geom_(geom) {
  oops::Log::trace() << "LocalUnstructuredInterpolator::LocalUnstructuredInterpolator start"
                     << std::endl;

  // note: creation of interpolators is posponed to apply / applyad calls
  // (because we don't know yet here which u/v/h masked/unmaked interpolators are required)

  oops::Log::trace() << "LocalUnstructuredInterpolator::LocalUnstructuredInterpolator done"
                     << std::endl;
  }

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::
apply(const oops::Variables & vars, const State & xx, const std::vector<bool> & mask,
      std::vector<double> & locvals) const {
  oops::Log::trace() << "LocalUnstructuredInterpolator::apply STATE start" << std::endl;

  atlas::FieldSet fset;
  xx.toFieldSet(fset, true);

  auto vals = locvals.begin();
  for (int i =0; i < vars.size(); i++) {
    auto interpolator = getInterpolator(vars[i]);

    // get a single variable
    oops::Variables var;
    var.push_back(vars[i]);

    // interpolate
    interpolator->apply(var, fset, mask, vals);
  }

  ASSERT(std::distance(locvals.begin(), vals) == locvals.size());
  oops::Log::trace() << "LocalUnstructuredInterpolator::apply STATE done" << std::endl;
}

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::
apply(const oops::Variables & vars, const Increment & dx, const std::vector<bool> & mask,
       std::vector<double> & locvals) const {
  oops::Log::trace() << "LocalUnstructuredInterpolator::apply Increment start" << std::endl;

  atlas::FieldSet fset;
  dx.toFieldSet(fset, true);

  auto vals = locvals.begin();
  for (int i =0; i < vars.size(); i++) {
    auto interpolator = getInterpolator(vars[i]);

    // get a single variable
    oops::Variables var;
    var.push_back(vars[i]);

    // interpolate
    interpolator->apply(var, fset, mask, vals);
  }

  oops::Log::trace() << "LocalUnstructuredInterpolator::apply Increment done" << std::endl;
}

// ------------------------------------------------------------------------------

void LocalUnstructuredInterpolator::
applyAD(const oops::Variables & vars, Increment & dx, const std::vector<bool> & mask,
        const std::vector<double> & locvals) const {
  oops::Log::trace() << "LocalUnstructuredInterpolator::applyAD start" << std::endl;

  int stride = locvals.size() / vars.size();

  Increment dz(dx);

  auto vals = locvals.begin();
  for (int i =0; i < vars.size(); i++) {
    auto interpolator = getInterpolator(vars[i]);

    // get a single variable
    oops::Variables var;
    var.push_back(vars[i]);

    dz.zero();
    atlas::FieldSet fset;
    dz.toFieldSet(fset, true);

    // interpolate
    interpolator->applyAD(var, fset, mask, vals);
    dx.getFieldSetAD(var, fset, false);
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

  // does the interpolator need to be created? (if it hasn't already yet)
  if (interp_[interp_idx].get() == nullptr) {
    // std::vector<double> lats_in;
    // std::vector<double> lons_in;
    // geom_.latlon(lats_in, lons_in, true, grid, masked);
    interp_[interp_idx] = std::make_shared<UnstructuredInterpolator>(
      config_, geom_, grid, masked, lats_out_, lons_out_);
  }

  // done, return the interpolator
  return interp_[interp_idx];
}

void LocalUnstructuredInterpolator::print(std::ostream & os) const {
}

}  // namespace soca
