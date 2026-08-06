/*
 * (C) Copyright 2024 UCAR.
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

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------

std::map<std::string, std::vector<std::string>> SocaLinVaderCookbook {
  {"sea_water_temperature", {"SeaWaterTemperature_A", "SeaWaterTemperature_B"}}
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

  default_ = (params_.linearVariableChangesWrapper.value() == boost::none);
}

// -----------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarTraj(const State & xfg,
                                         const oops::Variables & vars) {
  Log::trace() << "LinearVariableChange::changeVarTraj starting" << std::endl;

  // The following is TEMPORARY.
  // ----------------------------------------------------------------------------
  // We need to get some variables BEFORE we run VADER.
  State xfg_vadertraj(geom_, xfg);
  if (default_ && vars.has("sea_water_temperature")) {
    Log::debug() << "LinearVariableChange::changeVarTraj Pre-VADER variable changes. " << std::endl;
    oops::Variables preVaderVars(std::vector<std::string>{
      "latitude",
      "longitude",
      "sea_water_depth",
      "sea_area_fraction",
      });
    preVaderVars += xfg.variables();
    xfg_vadertraj.updateFields(preVaderVars);
    Model2GeoVaLs m2gv(geom_, params_.toConfiguration());
    m2gv.changeVar(xfg, xfg_vadertraj);
    Log::debug() << "LinearVariableChange::changeVarTraj variables after var change: "
        << xfg_vadertraj.variables() << std::endl;
  }

  // call Vader
  // ----------------------------------------------------------------------------
  oops::Variables varsFilled = xfg_vadertraj.variables();
  oops::Variables varsVader = vars;
  varsVader -= varsFilled;  // pass only the needed variables
  atlas::FieldSet xfs;
  xfg_vadertraj.toFieldSet(xfs);
  xfs.haloExchange();
  vader_->changeVarTraj(xfs, varsVader);

  // TODO(travis): this is not ideal. We are saving the first trajectory and
  // assuming it is the background. This should all get ripped out when the
  // variable changes that rely on the background are dealt with properly.
  if (!bkg_) { bkg_.reset(new State(xfg)); }

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
        LinearVariableChangeFactory::create(*bkg_, xfg, geom_, linVarChaPar));
    }
  } else {
    // No variable changes were specified, use the default (LinearModel2GeoVaLs)
    eckit::LocalConfiguration conf;
    conf.set("linear variable change name", "default");
    linVarChas_.push_back(LinearVariableChangeFactory::create(*bkg_, xfg, geom_,
      oops::validateAndDeserialize<GenericLinearVariableChangeParameters>(
        conf)));
  }
  Log::trace() << "LinearVariableChange::setTrajectory done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::initVaderTLAD(oops::Variables & ingredientVars) const {
  oops::Log::trace() << "LinearVariableChange::initVaderTLAD starting" << std::endl;
  oops::Variables originalIngredientVars = ingredientVars;
  varsVaderPopulates_ = vader_->initTLAD(ingredientVars);
  varsVaderPopulates_ -= originalIngredientVars;
  oops::Log::trace() << "LinearVariableChange::initVaderTLAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx,
                                       const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::changeVarTL starting" << std::endl;

  // If all variables already in incoming state just remove the no longer
  // needed fields
  // if (hasAllFields) {
  //   dx.updateFields(vars);
  //   oops::Log::trace() << "VariableChange::changeVar done (identity)"
  //                      << std::endl;
  //   return;
  // }

  // Make sure vader is fully initialized
  if (vader_->needsTLADInit()) {
    oops::Variables ingredientVars = dx.variables();
    initVaderTLAD(ingredientVars);
    oops::Log::info() << " changevarTL: vader will populate variables: "
                      << varsVaderPopulates_ << std::endl;
  }
  // If Vader is doing anything, call Vader
  if (varsVaderPopulates_.size() > 0) {
    atlas::FieldSet dxfs;
    dx.toFieldSet(dxfs);
    vader_->changeVarTL(dxfs);
    // Set intermediate state for the Increment containing original fields plus the ones
    // Vader has done
    oops::Variables varsVader = dx.variables();
    varsVader += varsVaderPopulates_;
    dx.updateFields(varsVader);
    dx.fromFieldSet(dxfs);
  }

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.validTime());

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
  Increment dxout(dx.geometry(), vars, dx.validTime());

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

  // Make sure vader is fully initialized
  if (vader_->needsTLADInit()) {
    oops::Variables ingredientVars(vars);
    initVaderTLAD(ingredientVars);
    oops::Log::info() << " changevarAD: vader will populate variables: "
                      << varsVaderPopulates_ << std::endl;
  }

  // Create dx_for_vader as a copy of dx, with just the variables created by Vader
  // (in the forward direction). Remove the variables created by vader from dx.
  // This way we ensure the model code will not be able to do the adjoint for these vars
  Increment dx_for_vader(dx, true);  // true => full copy
  oops::Variables varsVaderDidntPopulate = dx.variables();
  varsVaderDidntPopulate -= varsVaderPopulates_;
  dx_for_vader.updateFields(varsVaderPopulates_);
  dx.updateFields(varsVaderDidntPopulate);

  // Create empty output state
  Increment dxout(dx.geometry(), vars, dx.validTime());

  // Call model's adjoint variable change.
  for (ircst_ it = linVarChas_.rbegin(); it != linVarChas_.rend(); ++it) {
    dxout.zero();
    it->multiplyAD(dx, dxout);
    dx.updateFields(vars);
    dx = dxout;
  }

  // dxout needs to temporarily have the variables that Vader populated put into it before
  // being passed into vader_.changeVarAD, so Vader can do its adjoints.
  if (varsVaderPopulates_.size() > 0) {
    atlas::FieldSet dxfs;
    dx.toFieldSet(dxfs);
    oops::Variables varsVaderWillAdjoint = varsVaderPopulates_;
    atlas::FieldSet dxvaderfs;
    dx_for_vader.toFieldSet(dxvaderfs);
    for (const auto field : dxvaderfs) {
      dxfs.add(field);
    }

    oops::Variables varsAdjointed = vader_->changeVarAD(dxfs);
    varsVaderWillAdjoint -= varsAdjointed;
    // After changeVarAD, vader should have removed everything from varsVaderWillAdjoint,
    // indicating it did all the adjoints we expected it to.
    ASSERT(varsVaderWillAdjoint.size() == 0);
    // remove the fields that were adjointed by vader from dxout_fs
    util::removeFieldsFromFieldSet(dxfs, varsAdjointed.variables());
    // Copy dxout into dx for return
    dx.updateFields(vars);
    dx.fromFieldSet(dxfs);
  }

  Log::trace() << "LinearVariableChange::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx,
                                              const oops::Variables & vars) const {
  Log::trace() << "LinearVariableChange::multiplyInverseAD starting"
               << std::endl;

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.validTime());

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
