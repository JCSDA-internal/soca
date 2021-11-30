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

namespace soca {

// -------------------------------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom, const Parameters_ & params)
  : geom_(new Geometry(geom)), params_(params), linearVariableChange_(),
    fieldsMetadata_(geom.fieldsMetaData()) {}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::setTrajectory(const State & xbg, const State & xfg) {
  oops::Log::trace() << "LinearVariableChange::setTrajectory starting" << std::endl;
  // Create the variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(xbg, xfg, *geom_,
             params_.linearVariableChangeParametersWrapper.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiply(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiply starting" << std::endl;


  oops::Log::trace() << "LinearVariableChange::multiply done" << dx << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverse(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverse starting" << std::endl;


  oops::Log::trace() << "LinearVariableChange::multiplyInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyAD starting" << std::endl;


  oops::Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::multiplyInverseAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::multiplyInverseAD starting" << std::endl;


  oops::Log::trace() << "LinearVariableChange::multiplyInverseAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::print(std::ostream & os) const {
  os << "SOCA variable change";
}

// -------------------------------------------------------------------------------------------------

}  // namespace soca
