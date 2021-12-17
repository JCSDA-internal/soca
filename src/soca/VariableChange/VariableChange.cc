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
  Log::trace() << "VariableChange::changeVar starting" << vars << std::endl;

  // Create output state
  State xout(*x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVar(x, xout);

  x.updateFields(vars);
  // Copy data from temporary state
  x = xout;

// HOW CODE SHOULD LOOK  // Check whether vars already satisfied
// HOW CODE SHOULD LOOK  bool hasAllFields = x.hasAllFields();
// HOW CODE SHOULD LOOK
// HOW CODE SHOULD LOOK  if (hasAllFields) {
// HOW CODE SHOULD LOOK    x.updateFields(vars); // Remove any fields no longer needed
// HOW CODE SHOULD LOOK    Log::trace() << "VariableChange::changeVar done (identity)" << std::endl;
// HOW CODE SHOULD LOOK    return
// HOW CODE SHOULD LOOK  }
// HOW CODE SHOULD LOOK
// HOW CODE SHOULD LOOK  // Create output state
// HOW CODE SHOULD LOOK  State xout(*x.geometry(), vars, x.time());
// HOW CODE SHOULD LOOK
// HOW CODE SHOULD LOOK  // Call variable change
// HOW CODE SHOULD LOOK  variableChange_->changeVar(x, xout);
// HOW CODE SHOULD LOOK
// HOW CODE SHOULD LOOK  // Remove unused fields and allocate any new ones
// HOW CODE SHOULD LOOK  x.updateFields(vars);
// HOW CODE SHOULD LOOK
// HOW CODE SHOULD LOOK  // Copy data from temporary state
// HOW CODE SHOULD LOOK  x = xout;


  // Trace
  Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars)
  const {
  // Trace
  Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // Create output state
  State xout(*x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVarInverse(x, xout);

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
