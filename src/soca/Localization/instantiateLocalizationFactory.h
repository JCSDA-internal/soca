/*
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_LOCALIZATION_INSTANTIATELOCALIZATIONFACTORY_H_
#define SOCA_LOCALIZATION_INSTANTIATELOCALIZATIONFACTORY_H_

#include "soca/Traits.h"

#include "soca/Localization/Localization.h"

#include "oops/interface/LocalizationBase.h"

namespace soca {

void instantiateLocalizationFactory() {
}
// static oops::LocalizationMaker<soca::Traits,
//                                Localization> maker_("SOCA");
}  // namespace soca

#endif  // SOCA_LOCALIZATION_INSTANTIATELOCALIZATIONFACTORY_H_
