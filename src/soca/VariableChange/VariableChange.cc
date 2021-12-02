/*
 * (C) Copyright 2021 UCAR.
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
  // Trace
  Log::trace() << "VariableChange::changeVar starting" << std::endl;

  // Check if the incoming state has all the variables
  const bool hasAllFields = x.hasFields(vars);

  // If all variables already in incoming state just remove the no longer
  // needed fields
  // if (hasAllFields) {
  //   x.updateFields(vars);
  //   oops::Log::trace() << "VariableChange::changeVar done (identity)"
  //                      << std::endl;
  //   return;
  // }

  // Create output state
  State xout(*x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVar(x, xout);

  // Allocate any extra fields and remove fields no longer needed
  // x.updateFields(vars);

  // Copy data from temporary state
  x = xout;

  // Trace
  Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars)
  const {
  // Trace
  Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // Check if the incoming state has all the variables
  const bool hasAllFields = x.hasFields(vars);

  // If all variables already in incoming state just remove the no longer
  // needed fields
  // if (hasAllFields) {
  //   x.updateFields(vars);
  //   oops::Log::trace() << "VariableChange::changeVar done (identity)"
  //                      << std::endl;
  //   return;
  // }

  // Create output state
  State xout(*x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVarInverse(x, xout);

  // Allocate any extra fields and remove fields no longer needed
  // x.updateFields(vars);

  // Copy data from temporary state
  x = xout;

  // Trace
  Log::trace() << "VariableChange::changeVarInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::print(std::ostream & os) const {
  os << *variableChange_;
}

// -----------------------------------------------------------------------------

}  // namespace soca
