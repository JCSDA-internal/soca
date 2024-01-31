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

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// --------------------------------------------------------------------------------------
}  // namespace soca
