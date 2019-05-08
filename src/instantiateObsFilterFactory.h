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
#include "ufo/filters/BackgroundCheck.h"
#include "ufo/filters/BlackList.h"
#include "ufo/filters/ObsBoundsCheck.h"
#include "ufo/filters/ObsDomainCheck.h"
#include "ufo/filters/PreQC.h"
#include "ufo/filters/QCmanager.h"
#include "ufo/filters/Thinning.h"

namespace soca {

  void instantiateObsFilterFactory() {
    oops::instantiateObsFilterFactory<Traits>();
    static oops::FilterMaker<Traits,
      oops::ObsFilter<Traits, ufo::QCmanager>>
      makerChkManager_("QCmanager");
    static oops::FilterMaker<Traits,
      oops::ObsFilter<Traits, ufo::PreQC>>
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
    static oops::FilterMaker<Traits,
      oops::ObsFilter<Traits, ufo::Thinning>>
      makerThinning_("Thinning");
  }
}  // namespace soca

#endif  // SRC_INSTANTIATEOBSFILTERFACTORY_H_
