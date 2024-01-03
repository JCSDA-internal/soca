/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <map>
#include <ostream>
#include <string>

#include "soca/LinearVariableChange/LinearVariableChange.h"
#include "soca/VariableChange/Model2GeoVaLs/Model2GeoVaLs.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------

std::map<std::string, std::vector<std::string>> SocaLinVaderCookbook {
  {"sea_water_temperature", {"SeaWaterTemperature_A"}},
};

// -----------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom,
                                           const eckit::Configuration & config)
  : geom_(geom), params_(), linVarChas_(), vader_(), default_(false) {
    params_.deserialize(config);

  // setup vader
  eckit::LocalConfiguration vaderConfig, vaderCookbookConfig;
  for (auto kv : SocaLinVaderCookbook) vaderCookbookConfig.set(kv.first, kv.second);
  vaderConfig.set(vader::configCookbookKey, vaderCookbookConfig);
  vader_.reset(new vader::Vader(params_.vader, vaderConfig));

  const boost::optional<std::vector<LinearVariableChangeParametersWrapper>> &
    linVarChgs = params_.linearVariableChangesWrapper;
  default_ = linVarChgs == boost::none;
}

// -----------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarTraj(const State & xfg,
                                         const oops::Variables & vars) {
  Log::trace() << "LinearVariableChange::changeVarTraj starting, default: "<<default_ << std::endl;

  Log::debug() << "LinearVariableChange::changeVarTraj vars in: "
               << xfg.variables() << std::endl;
  Log::debug() << "LinearVariableChange::changeVarTraj vars out: "
               << vars << std::endl;
  
  // The following is TEMPORARY.
  // ----------------------------------------------------------------------------
  // We need to do some variable renaming BEFORE we run VADER.
  // Eventually, we will internally rename these variables when they are
  // first loaded in so that we won't have to worry about it here.
  State xfg2(xfg);
  if (default_ && vars.has("sea_water_temperature")) {
    Log::debug() << "LinearVariableChange::changeVarTraj Pre-VADER variable changes. " << std::endl;
    oops::Variables preVaderVars(std::vector<std::string>{
      "latitude",
      "longitude",
      "sea_water_potential_temperature",
      "sea_water_salinity",
      "sea_water_depth",
      "sea_area_fraction",
      });
    preVaderVars += xfg.variables();
    xfg2.updateFields(preVaderVars);
    Model2GeoVaLs m2gv(xfg.geometry(), params_.toConfiguration());
    m2gv.changeVar(xfg, xfg2);
    Log::debug() << "LinearVariableChange::changeVarTraj variables after var change: "
        << xfg2.variables() << std::endl;
  }

  // call Vader
  // ----------------------------------------------------------------------------
  Log::debug() << "LinearVariableChange::changeVarTraj VADER variable changes. " << std::endl;
  oops::Variables varsFilled = xfg2.variables();
  oops::Variables varsVader = vars;
  varsVader -= varsFilled; // pass only the needed variables
  atlas::FieldSet xfs;
  xfg2.toFieldSet(xfs);
  varsVaderPopulates_ = vader_->changeVarTraj(xfs, varsVader);
  varsFilled += varsVaderPopulates_;
  xfg2.updateFields(varsFilled);
  xfg2.fromFieldSet(xfs);
  Log::debug() << "LinearVariableChange::changeVarTraj variables after var change: "
    << xfg2.variables() << std::endl;

  // TODO(travis): this is not ideal. We are saving the first trajectory and
  // assuming it is the background. This should all get ripped out when the
  // variable changes that rely on the background are dealt with properly.
  if (!bkg_) { bkg_.reset(new State(xfg2)); }

  const boost::optional<std::vector<LinearVariableChangeParametersWrapper>> &
    linVarChgs = params_.linearVariableChangesWrapper;

  if (linVarChgs != boost::none) {
    // If one or more linear variable changes were specified:

    // Create the linear variable change(s)
    for (const LinearVariableChangeParametersWrapper & linVarChaParWra :
        *linVarChgs) {
      // Get parameters for this linear variable change
      const LinearVariableChangeParametersBase & linVarChaPar =
            linVarChaParWra.linearVariableChangeParameters;
      // Add linear variable change to vector
      linVarChas_.push_back(
        LinearVariableChangeFactory::create(*bkg_, xfg2, geom_, linVarChaPar));
    }
  } else {  
    // No variable changes were specified, use the default (LinearModel2GeoVaLs)
    eckit::LocalConfiguration conf;
    conf.set("linear variable change name", "default");
    linVarChas_.push_back(LinearVariableChangeFactory::create(*bkg_, xfg2, geom_,
      oops::validateAndDeserialize<GenericLinearVariableChangeParameters>(
        conf)));
  }
  Log::trace() << "LinearVariableChange::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx,
                                       const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::changeVarTL starting, default: " << default_ << std::endl;

  Log::debug() << "LinearVariableChange::changeVarTL vars in: "
               << dx.variables() << std::endl;
  Log::debug() << "LinearVariableChange::changeVarTL vars out: "
               << vars << std::endl;


  // The following is TEMPORARY.
  // ----------------------------------------------------------------------------
  // We need to do some variable renaming BEFORE we run VADER.
  // Eventually, we will internally rename these variables when they are
  // first loaded in so that we won't have to worry about it here.
  if (default_ && vars.has("sea_water_temperature")) {
    Log::debug() << "LinearVariableChange::changeVarTL Pre-VADER variable changes. " << std::endl;
    oops::Variables preVaderVars(std::vector<std::string>{
      "latitude",
      "longitude",
      "sea_water_potential_temperature",
      "sea_water_salinity",
      "sea_water_depth",
      "sea_area_fraction"});
    preVaderVars += dx.variables();    
    Increment preVader(dx.geometry(), preVaderVars, dx.time());

    for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
      it->multiply(dx, preVader);
      dx.updateFields(preVaderVars);
      dx = preVader;
    }
    Log::debug() << "LinearVariableChange::changeVarTL variables after var change: "
                << dx.variables() << std::endl;  
  }


  // call Vader
  // ----------------------------------------------------------------------------
  Log::debug() << "LinearVariableChange::changeVarTL VADER variable changes. " << std::endl;
  atlas::FieldSet xfs;
  dx.toFieldSet(xfs);
  oops::Variables varsFilled = dx.variables();
  oops::Variables varsVader = vars;
  varsVader -= varsFilled;
  varsFilled += vader_->changeVarTL(xfs, varsVader);
  dx.updateFields(varsFilled);
  dx.fromFieldSet(xfs);
  Log::debug() << "LinearVariableChange::changeVarTL variables after var change: "
               << dx.variables() << std::endl;

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
     dxout.zero();
     it->multiply(dx, dxout);
     dx.updateFields(vars);
     dx = dxout;
  }

  Log::trace() << "LinearVariableChange::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseTL(Increment & dx,
                                              const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyInverse starting"
               << vars << std::endl;

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (ircst_ it = linVarChas_.rbegin(); it != linVarChas_.rend(); ++it) {
    dxout.zero();
    it->multiplyInverse(dx, dxout);
    dx.updateFields(vars);
    dx = dxout;
  }

  Log::trace() << "LinearVariableChange::multiplyInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarAD(Increment & dx,
                                       const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::changeVarAD starting" << std::endl;

  Log::debug() << "LinearVariableChange::changeVarAD vars in: "
               << dx.variables() << std::endl;
  Log::debug() << "LinearVariableChange::changeVarAD vars out: "
               << vars << std::endl;

  // NOTE: the IF is temporary, we need to do some variable renaming afterward
  if (default_ && dx.variables().has("sea_water_temperature")) {    
    // call vader    
    Log::debug() << "LinearVariableChange::changeVarAD VADER variable changes. " << std::endl;    
    oops::Variables varsToAdj(std::vector<std::string>{
      "sea_water_potential_temperature"});
    oops::Variables varsToDrop(std::vector<std::string>{
      "sea_water_temperature"});
    
    // run Vader    
    oops::Variables vaderVars = dx.variables();
    vaderVars += varsToAdj;
    dx.updateFields(vaderVars);
  
    atlas::FieldSet xfs;
    dx.toFieldSet(xfs);

    oops::Variables varsVaderWillAdjoint = varsVaderPopulates_;
    vader_->changeVarAD(xfs, varsVaderWillAdjoint);
    ASSERT(varsVaderWillAdjoint.size() == 0);

    vaderVars -= varsToDrop;
    dx.updateFields(vaderVars);
    dx.fromFieldSet(xfs);
    Log::debug() << "LinearVariableChange::changeVarTL variables after var change: "
                 << dx.variables() << std::endl;  
  }


  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (ircst_ it = linVarChas_.rbegin(); it != linVarChas_.rend(); ++it) {
    dxout.zero();
    it->multiplyAD(dx, dxout);
    dx.updateFields(vars);
    dx = dxout;
  }

  Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx,
                                              const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyInverseAD starting"
               << std::endl;

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change(s)
  for (icst_ it = linVarChas_.begin(); it != linVarChas_.end(); ++it) {
    dxout.zero();
    it->multiplyInverseAD(dx, dxout);
    dx.updateFields(vars);
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
