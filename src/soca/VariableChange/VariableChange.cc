/*
 * (C) Copyright 2021-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"
#include "soca/VariableChange/VariableChange.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------

VariableChange::VariableChange(const Parameters_ & params,
                               const Geometry & geometry) {
  // Create the variable change
  variableChange_.reset(VariableChangeFactory::create(geometry,
      params.variableChangeParametersWrapper.variableChangeParameters.value()));
}

// -----------------------------------------------------------------------------

VariableChange::~VariableChange() {}

// -----------------------------------------------------------------------------

void VariableChange::changeVar(State & x, const oops::Variables & vars) const {
  Log::trace() << "VariableChange::changeVar starting" << std::endl;

  Log::debug() << "VariableChange::changeVar vars in: "
               << x.variables() << std::endl;
  Log::debug() << "VariableChange::changeVar vars out: "
               << vars << std::endl;

  // TODO(travis) rename in/out variables so that skipping this
  // works for Model2Ana (i.e. we need rotated/unrotate u/v renamed different)
  // // If the variables are the same, don't bother doing anything!
  // if (!(x.variables() == vars)) {

  // Create output state
  State xout(x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVar(x, xout);

  // Copy data from temporary state
  x.updateFields(vars);
  x = xout;

  // }

  Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x,
                                      const oops::Variables & vars) const {
  Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  Log::debug() << "VariableChange::changeVarInverse vars in: "
               << x.variables() << std::endl;
  Log::debug() << "VariableChange::changeVarInverse vars out: "
               << vars << std::endl;

  // If the variables are the same, don't bother doing anything!
  if (!(x.variables() == vars)) {
    // Create output state
    State xout(x.geometry(), vars, x.time());

    // Call variable change
    variableChange_->changeVarInverse(x, xout);

    // Copy data from temporary state
    x.updateFields(vars);
    x = xout;
  }

  Log::trace() << "VariableChange::changeVarInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::print(std::ostream & os) const {
  os << *variableChange_;
}

// -----------------------------------------------------------------------------

}  // namespace soca
