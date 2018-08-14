/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_INSTANTIATECHANGEVARFACTORY_H_
#define SOCA_MODEL_INSTANTIATECHANGEVARFACTORY_H_

#include "src/BalanceOperators/Kst/Kst.h"
#include "src/Traits.h"
#include "oops/interface/LinearVariableChange.h"

namespace soca {

void instantiateBalanceOpFactory() {
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::Kst> >
              makerBalanceOpSOCA_("KstSOCA");
  /*static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::Ksshts> >
              makerBalanceOpSOCA_("KsshtsSOCA");*/
}

}  // namespace qg

#endif  // SOCA_MODEL_INSTANTIATECHANGEVARFACTORY_H_
