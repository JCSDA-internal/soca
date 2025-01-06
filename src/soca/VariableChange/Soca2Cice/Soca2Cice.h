/*
 * (C) Copyright 2022-2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "soca/Traits.h"

#include "soca/VariableChange/Base/VariableChangeBase.h"

namespace soca {

class Soca2Cice: public VariableChangeBase {
 public:
  // ---------------------------------------------------------------------------
  /// @brief Parameters for rescaling CICE analysis in the ice pack
  class RescaleParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(RescaleParameters, oops::Parameters)
   public:
    oops::Parameter<bool> rescale{"rescale",
      "rescale analysis in the ice pack", false, this};
    oops::Parameter<float> min_hice{"min hice",
      "min ice thickness to trigger adjusting ice volume", 0.5, this};
    oops::Parameter<float> min_hsno{"min hsno",
      "min snow thickness to trigger adjusting snow volume", 0.1, this};
  };

  /// @brief Parameters for the redistribution of sea ice
  class RedistributionParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(RedistributionParameters, oops::Parameters)
   public:
    oops::Parameter<double> edge{"seaice edge",
      "Threshold for sea ice edge, used in shuffle and rescale prior options",
      0.15, this, {oops::minConstraint(0.0), oops::maxConstraint(1.0)}};
    oops::Parameter<bool> shuffle{"shuffle",
      "Option to shuffle sea ice in the marginal ice zone (where ice concentration"
      " < seaice edge)",
      true, this};
    oops::Parameter<RescaleParameters> rescale{"rescale prior",
      "Option to rescale sea ice in the ice pack zone (where ice concentration"
      " > seaice edge)", {}, this};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for CICE restart background
  class BackgroundParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(BackgroundParameters, oops::Parameters)
   public:
    oops::RequiredParameter<std::string> restart{"restart", this};
    oops::RequiredParameter<int> ncat{"ncat", this, {oops::minConstraint(1)}};
    oops::RequiredParameter<int> ice_lev{"ice_lev", this, {oops::minConstraint(1)}};
    oops::RequiredParameter<int> sno_lev{"sno_lev", this, {oops::minConstraint(1)}};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for CICE restart output
  class OutputParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(OutputParameters, oops::Parameters)
   public:
    oops::RequiredParameter<std::string> restart{"restart", this};
  };

  // ---------------------------------------------------------------------------
  /// @brief Parameters for adding soca increment to CICE restart files
  class Parameters : public VariableChangeParametersBase {
    OOPS_CONCRETE_PARAMETERS(Parameters, VariableChangeParametersBase)
   public:
    oops::RequiredParameter<BackgroundParameters> background{"cice background state", this};
    oops::RequiredParameter<OutputParameters> output{"cice output", this};
    oops::OptionalParameter<eckit::LocalConfiguration> incOutput{"increment output", this};
    oops::OptionalParameter<eckit::LocalConfiguration> incInput{"soca increment", this};
    oops::RequiredParameter<RedistributionParameters> arctic{"arctic", this};
    oops::RequiredParameter<RedistributionParameters> antarctic{"antarctic", this};
  };

  const std::string classname() {return "soca::Soca2Cice";}

  Soca2Cice(const Geometry &, const eckit::Configuration &);
  ~Soca2Cice();

  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  const Geometry & geom_;
  Parameters params_;
  int keySoca2Cice_;
  void print(std::ostream &) const override {}
};

}  // namespace soca
