/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "oops/base/ModelBase.h"
#include "oops/interface/Model.h"

#include "soca/Model/OceanIceEmulator/ModelOceanIceEmulator.h"
#include "soca/Traits.h"

namespace soca {

void Traits::registerModels() {
  static oops::ModelMaker<Traits, oops::interface::Model<Traits, ModelOceanIceEmulator> >
      makerSoca_("ModelOceanIceEmulator");
}

}  // namespace soca
