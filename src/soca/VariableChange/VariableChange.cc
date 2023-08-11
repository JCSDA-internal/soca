/*
 * (C) Copyright 2021-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>
#include <vector>

#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"
#include "soca/VariableChange/VariableChange.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------

VariableChange::VariableChange(const Parameters_ & params,
                               const Geometry & geometry) {
  // setup vader
  vader_.reset(new vader::Vader(params.vader));

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
  // if (!(x.variables() == vars))

  // The following is TEMPORARY.
  // ----------------------------------------------------------------------------
  // We need to do some variable renaming BEFORE we run VADER.
  // Eventually, we will internally rename these variables when they are
  // first loaded in so that we won't have to worry about it here.
  if (vars.has("sea_water_temperature")) {
    Log::debug() << "VariableChange::changeVar Pre-VADER variable changes. " << std::endl;
    oops::Variables preVaderVars(std::vector<std::string>{
      "latitude",
      "sea_water_potential_temperature",
      "sea_water_salinity",
      "sea_water_depth"});
    preVaderVars += x.variables();
    State preVader(x.geometry(), preVaderVars, x.time());
    variableChange_->changeVar(x, preVader);
    x.updateFields(preVaderVars);
    x = preVader;
    Log::debug() << "VariableChange::changeVar variables after var change: "
                << x.variables() << std::endl;
  }

  // call Vader
  // ----------------------------------------------------------------------------
  Log::debug() << "VariableChange::changeVar VADER variable changes. " << std::endl;
  // Record start variables
  oops::Variables varsFilled = x.variables();
  oops::Variables varsVader = vars;
  varsVader -= varsFilled;  // Pass only the needed variables

  // Call Vader. On entry, varsVader holds the vars requested from Vader; on exit,
  // it holds the vars NOT fulfilled by Vader, i.e., the vars still to be requested elsewhere.
  // vader_->changeVar also returns the variables fulfilled by Vader. These variables are allocated
  // and populated and added to the FieldSet (xfs).
  atlas::FieldSet xfs;
  x.toFieldSet(xfs);
  varsFilled += vader_->changeVar(xfs, varsVader);
  x.updateFields(varsFilled);
  x.fromFieldSet(xfs);
  Log::debug() << "VariableChange::changeVar variables after var change: "
               << x.variables() << std::endl;

  // soca specific transforms
  // ----------------------------------------------------------------------------
  Log::debug() << "VariableChange::changeVar SOCA specific post-VADER variable changes. "
               << std::endl;
  State xout(x.geometry(), vars, x.time());
  variableChange_->changeVar(x, xout);
  x.updateFields(vars);
  x = xout;
  Log::debug() << "VariableChange::changeVar variables after var change: "
               << x.variables() << std::endl;


  Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -----------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x,
                                      const oops::Variables & vars) const {
  changeVar(x, vars);
}

// -----------------------------------------------------------------------------

void VariableChange::print(std::ostream & os) const {
  os << *variableChange_;
}

// -----------------------------------------------------------------------------

}  // namespace soca
