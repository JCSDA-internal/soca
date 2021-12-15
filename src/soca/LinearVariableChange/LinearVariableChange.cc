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

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    dxout.zero();
    it->multiply(dx, dxout);
    dx.updateFields(dxout);
  }

  Log::trace() << "LinearVariableChange::multiply done" << dx << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::multiplyInverse(Increment & dx,
                                           const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyInverse starting" << std::endl;

  // BUG vars points to nowhere ...
  return;
  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    std::cout << " --------------- dx:" << dx << std::endl;
    dxout.zero();
    it->multiplyInverse(dx, dxout);
    dx.updateFields(dxout);
  }

  Log::trace() << "LinearVariableChange::multiplyInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::multiplyAD(Increment & dx,
                                           const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyAD starting" << std::endl;
  std::cout << "------ vars in linvarchange multiplyad: " << vars << std::endl;

  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    dxout.zero();
    it->multiplyAD(dx, dxout);
    dx.updateFields(dxout);
  }

  Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::multiplyInverseAD(Increment & dx,
                                          const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyInverseAD starting"
               << std::endl;
  return;
  // Create output state
  Increment dxout(*dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    dxout.zero();
    it->multiplyInverseAD(dx, dxout);
    dx.updateFields(dxout);
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
