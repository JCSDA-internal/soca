/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_INSTANTIATEBALANCEOPFACTORY_H_
#define SOCA_TRANSFORMS_INSTANTIATEBALANCEOPFACTORY_H_

#include "soca/Transforms/Balance/Balance.h"
#include "soca/Transforms/BkgErr/BkgErr.h"
#include "soca/Transforms/BkgErrGodas/BkgErrGodas.h"
#include "soca/Transforms/BkgErrFilt/BkgErrFilt.h"
#include "soca/Transforms/HorizFilt/HorizFilt.h"
#include "soca/Transforms/VertConv/VertConv.h"
#include "soca/Traits.h"
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
  static oops::LinearVariableChangeMaker<soca::Traits,
              oops::LinearVariableChange<soca::Traits, soca::HorizFilt> >
              makerBalanceOpHorizFILT_("HorizFiltSOCA");
}
}  // namespace soca

#endif  // SOCA_TRANSFORMS_INSTANTIATEBALANCEOPFACTORY_H_
