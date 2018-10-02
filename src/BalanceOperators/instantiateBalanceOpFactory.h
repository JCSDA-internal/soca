/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_INSTANTIATECHANGEVARFACTORY_H_
#define SOCA_MODEL_INSTANTIATECHANGEVARFACTORY_H_

#include "src/BalanceOperators/Kst/Kst.h"
#include "src/BalanceOperators/Ktc/Ktc.h"
#include "src/BalanceOperators/BkgErr/BkgErr.h"
#include "src/BalanceOperators/Ksshts/Ksshts.h"
#include "src/BalanceOperators/VertConv/VertConv.h"
#include "src/Traits.h"
#include "oops/interface/LinearVariableChange.h"

namespace soca {

void instantiateBalanceOpFactory() {
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::VertConv> >
              makerBalanceOpVertConvSOCA_("VertConvSOCA");  
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::Kst> >
              makerBalanceOpKstSOCA_("KstSOCA");
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::Ksshts> >
              makerBalanceOpKsshtsSOCA_("KsshtsSOCA");
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::Ktc> >
              makerBalanceOpKtcSOCA_("KtcSOCA");
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::BkgErr> >
              makerBalanceOpBkgErrSOCA_("BkgErrSOCA");  
}

}  // namespace qg

#endif  // SOCA_MODEL_INSTANTIATECHANGEVARFACTORY_H_
