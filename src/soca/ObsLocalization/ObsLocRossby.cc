/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/ObsLocalization/ObsLocRossby.h"
#include "soca/Traits.h"

#include "oops/base/ObsLocalizationBase.h"
#include "ufo/ObsTraits.h"


namespace soca {

static oops::ObsLocalizationMaker< soca::Traits, ufo::ObsTraits,
    soca::ObsLocRossby<soca::Traits>> obslocrossby_("Rossby");

}  // namespace soca
