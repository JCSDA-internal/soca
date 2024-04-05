/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "saber/blocks/SaberBlockParametersBase.h"

namespace soca {

// --------------------------------------------------------------------------------------

class ExplicitDiffusionParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionParameters, saber::SaberBlockParametersBase)
 public:
  oops::OptionalParameter<eckit::LocalConfiguration> read{"read", this};
  oops::OptionalParameter<eckit::LocalConfiguration> calibration{"calibration", this};

  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};
  oops::OptionalParameter<eckit::LocalConfiguration> groups{"group mapping", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// --------------------------------------------------------------------------------------
}  // namespace soca
