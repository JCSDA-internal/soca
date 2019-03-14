/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_INSTANTIATEOBSFILTERFACTORY_H_
#define SRC_INSTANTIATEOBSFILTERFACTORY_H_

#include "src/Traits.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/BackgroundCheck.h"
#include "ufo/BlackList.h"
#include "ufo/ObsBoundsCheck.h"
#include "ufo/ObsDomainCheck.h"
#include "ufo/ObsPreQC.h"

namespace soca {

  void instantiateObsFilterFactory() {
    oops::instantiateObsFilterFactory<Traits>();
    static oops::FilterMaker<Traits, oops::ObsFilter<Traits, ufo::ObsPreQC>>
      makerChkpreqc_("PreQC");
    static oops::FilterMaker<Traits,
      oops::ObsFilter<Traits, ufo::ObsDomainCheck>>
      makerChkdomaincheck_("Domain Check");
    static oops::FilterMaker<Traits,
      oops::ObsFilter<Traits, ufo::ObsBoundsCheck>>
      makerChkbndcheck_("Bounds Check");
    static oops::FilterMaker<Traits,
      oops::ObsFilter<Traits, ufo::BlackList>>
      makerChkblcklist_("BlackList");
    static oops::FilterMaker<Traits,
      oops::ObsFilter<Traits, ufo::BackgroundCheck>>
      makerBkgChk_("Background Check");
  }
}  // namespace soca

#endif  // SRC_INSTANTIATEOBSFILTERFACTORY_H_
