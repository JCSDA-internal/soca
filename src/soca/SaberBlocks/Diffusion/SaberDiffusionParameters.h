/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "saber/blocks/SaberBlockParametersBase.h"

namespace soca {
  
class SaberDiffusionParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(SaberDiffusionParameters, saber::SaberBlockParametersBase)
 public:
  oops::OptionalParameter<eckit::LocalConfiguration> read{"read", this};
  oops::OptionalParameter<eckit::LocalConfiguration> calibration{"calibration", this};
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

}