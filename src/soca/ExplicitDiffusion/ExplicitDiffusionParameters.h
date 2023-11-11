/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "saber/blocks/SaberBlockParametersBase.h"

namespace soca {

// --------------------------------------------------------------------------------------

class ExplicitDiffusionIOParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionIOParameters, oops::Parameters)
 public:
  oops::RequiredParameter<std::string> filename{"filename", this};
  oops::OptionalParameter<std::string> varName{"variable name", this};
};

// --------------------------------------------------------------------------------------

class ExplicitDiffusionScalesParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionScalesParameters, oops::Parameters)
 public:
  oops::OptionalParameter<double> fixedValue{"fixed value", this};
  oops::OptionalParameter<ExplicitDiffusionIOParameters> fromFile{"from file", this};
  oops::Parameter<bool> asGaussian{"as gaussian", false, this};
};

// --------------------------------------------------------------------------------------

class ExplicitDiffusionCalibrationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionCalibrationParameters, oops::Parameters)
 public:
  oops::RequiredParameter<eckit::LocalConfiguration> normalization{"normalization", this};
  oops::RequiredParameter<ExplicitDiffusionScalesParameters> scales{"scales", this};
  oops::OptionalParameter<ExplicitDiffusionIOParameters> write {"write", this};
};

// --------------------------------------------------------------------------------------

class ExplicitDiffusionParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionParameters, saber::SaberBlockParametersBase)
 public:
  oops::OptionalParameter<ExplicitDiffusionIOParameters> read{"read", this};
  oops::OptionalParameter<ExplicitDiffusionCalibrationParameters>
    calibration{"calibration", this};
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// --------------------------------------------------------------------------------------
}  // namespace soca
