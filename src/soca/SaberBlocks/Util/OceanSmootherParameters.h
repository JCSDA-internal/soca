/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace soca
{

class OceanSmootherParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OceanSmootherParameters, oops::Parameters)

 public:

  // ----------------------------------------------------------------------------------------------
  class Horizontal : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Horizontal, oops::Parameters)
   public:
    oops::Parameter<std::string> rossbyVariable{"rossby radius variable", "rossby_radius", this};

    oops::Parameter<double> base{"base value", 0.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> rossbyMult{"rossby mult", 1.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> minGridMult{"min grid mult", 1.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> min{"min", 0.0, this, {oops::minConstraint(0.0)}};
    oops::Parameter<double> max{"max", std::numeric_limits<double>::max(), this,
                                {oops::minConstraint(0.0)}};
  };

  // ----------------------------------------------------------------------------------------------

  class Vertical : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(Vertical, oops::Parameters)
   public:
    oops::RequiredParameter<size_t> levels{"levels", this};
    oops::Parameter<double> base{"base value", 1.0, this, {oops::minConstraint(0.0)}};
  };

  // ----------------------------------------------------------------------------------------------

  oops::Parameter<std::string> mask{"mask variable", "", this};
  oops::OptionalParameter<Horizontal> horizontal{"horizontal", this};
  oops::OptionalParameter<Vertical> vertical{"vertical", this};
};

}  // namespace soca
