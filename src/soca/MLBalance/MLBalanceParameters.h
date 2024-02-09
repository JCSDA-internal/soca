#pragma once

#include <string>
#include <vector>

#include "saber/blocks/SaberBlockParametersBase.h"

namespace soca {

// --------------------------------------------------------------------------------------

class MLBalanceParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MLBalanceParameters, saber::SaberBlockParametersBase)
 public:

  oops::RequiredParameter<eckit::LocalConfiguration> mlbalances{"ML Balances", this};
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};
  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables({"tocn", "socn", "ssh", "cicen", "hicen", "hsnon"});}
};

// --------------------------------------------------------------------------------------
}  // namespace soca
