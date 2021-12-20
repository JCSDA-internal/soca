/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_OBSLOCALIZATION_OBSLOCROSSBYPARAMETERS_H_
#define SOCA_OBSLOCALIZATION_OBSLOCROSSBYPARAMETERS_H_

#include "ufo/obslocalization/ObsHorLocParameters.h"

namespace soca {

class ObsLocRossbyParameters : public ufo::ObsHorLocParameters {
  OOPS_CONCRETE_PARAMETERS(ObsLocRossbyParameters, ufo::ObsHorLocParameters)

 public:
  oops::Parameter<double> base{"base value", 0.0, this};
  oops::Parameter<double> mult{"rossby mult", 1.0, this};
  oops::Parameter<double> min_grid{"min grid mult", 1.0, this};
  oops::OptionalParameter<double> min{"min value", this};
  oops::OptionalParameter<double> max{"max value", this};
};

}  // namespace soca
#endif  // SOCA_OBSLOCALIZATION_OBSLOCROSSBYPARAMETERS_H_
