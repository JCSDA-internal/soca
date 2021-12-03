/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "soca/LinearVariableChange/LinearVariableChange.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom,
                                           const Parameters_ & params)
  : geom_(new Geometry(geom)), params_(params), linVarChas_() {}

// -----------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -----------------------------------------------------------------------------

void LinearVariableChange::setTrajectory(const State & xbg, const State & xfg) {
  Log::trace() << "LinearVariableChange::setTrajectory starting" << std::endl;
  // Create the linear variable change(s)
  for (const LinearVariableChangeParametersWrapper & linVarChaParWra :
       params_.linearVariableChangesWrapper.value()) {
    // Get parameters for this linear variable change
    const LinearVariableChangeParametersBase & linVarChaPar =
           linVarChaParWra.linearVariableChangeParameters;
    // Add linear variable change to vector
    linVarChas_.push_back(LinearVariableChangeFactory::create(xbg, xfg, *geom_,
                          linVarChaPar));
  }
  Log::trace() << "LinearVariableChange::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::multiply(Increment & dx,
                                    const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiply starting" << std::endl;

  // Check if the incoming state has all the variables
  const bool hasAllFields = dx.hasFields(vars);

  // If all variables already in incoming state just remove the no longer
  // needed fields
  // if (hasAllFields) {
  //   dx.updateFields(vars);
  //   oops::Log::trace() << "VariableChange::changeVar done (identity)"
  //                      << std::endl;
  //   return;
  // }

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    dxout.zero();
    it->multiply(dx, dxout);
    dx = dxout;
  }

  // Allocate any extra fields and remove fields no longer needed
  // dx.updateFields(vars);

  // Copy data from temporary state
  // dx = dxout;

  Log::trace() << "LinearVariableChange::multiply done" << dx << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::multiplyInverse(Increment & dx,
                                           const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyInverse starting" << vars << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (ircst_ it = linVarChas_.rbegin(); it != linVarChas_.rend(); ++it) {
    dxout.zero();
    it->multiplyInverse(dx, dxout);
    dx = dxout;
  }

  Log::trace() << "LinearVariableChange::multiplyInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::multiplyAD(Increment & dx,
                                           const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyAD starting" << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (ircst_ it = linVarChas_.rbegin(); it != linVarChas_.rend(); ++it) {
    dxout.zero();
    it->multiplyAD(dx, dxout);
    dx = dxout;
  }

  Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::multiplyInverseAD(Increment & dx,
                                          const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyInverseAD starting"
               << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    dxout.zero();
    it->multiplyInverseAD(dx, dxout);
    dx = dxout;
  }

  Log::trace() << "LinearVariableChange::multiplyInverseAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::print(std::ostream & os) const {
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    os << *it;
  }
}

// -----------------------------------------------------------------------------

}  // namespace soca
