/*
 * (C) Copyright 2025-2025 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include "oops/base/LinearModelBase.h"
#include "oops/base/ModelBase.h"
#include "oops/interface/LinearModelInterface.h"
#include "oops/interface/ModelInterface.h"

#include "soca/LinearModel/OceanIceEmulator/LinearModelOceanIceEmulator.h"
#include "soca/Model/OceanIceEmulator/ModelOceanIceEmulator.h"

namespace soca {

template <typename MODEL>
void instantiateModelFactory() {
  static oops::ModelMaker<MODEL, oops::ModelInterface<MODEL, ModelOceanIceEmulator> >
                    makerModelOIE_("ModelOceanIceEmulator");

  static oops::LinearModelMaker<MODEL,
                                oops::LinearModelInterface<MODEL, LinearModelOceanIceEmulator> >
                     makerTLMOIE_("LinearModelOceanIceEmulator");
}

}  // namespace soca

