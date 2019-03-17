/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_TRANSFORMS_INSTANTIATEBALANCEOPFACTORY_H_
#define SRC_TRANSFORMS_INSTANTIATEBALANCEOPFACTORY_H_

#include "src/Transforms/Balance/Balance.h"
#include "src/Transforms/BkgErr/BkgErr.h"
#include "src/Transforms/BkgErrGodas/BkgErrGodas.h"
#include "src/Transforms/BkgErrFilt/BkgErrFilt.h"
#include "src/Transforms/VertConv/VertConv.h"
#include "src/Traits.h"
#include "oops/interface/LinearVariableChange.h"

namespace soca {

void instantiateBalanceOpFactory() {
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::VertConv> >
              makerBalanceOpVertConvSOCA_("VertConvSOCA");
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::BkgErr> >
              makerBalanceOpBkgErrSOCA_("BkgErrSOCA");
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::BkgErrGodas> >
              makerBalanceOpBkgErrGODAS_("BkgErrGODAS");
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::BkgErrFilt> >
              makerBalanceOpBkgErrFILT_("BkgErrFILT");
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::Balance> >
              makerBalanceOpBalanceSOCA_("BalanceSOCA");
}
}  // namespace soca

#endif  // SRC_TRANSFORMS_INSTANTIATEBALANCEOPFACTORY_H_
