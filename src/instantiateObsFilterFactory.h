/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_INSTANTIATEOBSFILTERFACTORY_H_
#define SOCA_MODEL_INSTANTIATEOBSFILTERFACTORY_H_

#include "src/Traits.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/BackgroundCheck.h"
#include "ufo/qcfilters/domaincheck/DomainCheck.h"

namespace soca {

  void instantiateObsFilterFactory() {
    oops::instantiateObsFilterFactory<Traits>();
    static oops::FilterMaker<Traits,oops::ObsFilter<Traits, ufo::BackgroundCheck>>
      makerBkgChk_("Background Check");
    static oops::FilterMaker<Traits,oops::ObsFilter<Traits, ufo::DomainCheck>>
      makerDmChk_("Domain Check");
  }

}  // namespace soca

#endif  // SOCA_MODEL_INSTANTIATEOBSFILTERFACTORY_H_
