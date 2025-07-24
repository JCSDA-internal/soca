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

#include "mpasjedi/Model/Model.h"
#include "mpasjedi/Tlm/Tlm.h"

namespace soca {

template <typename MODEL>
void instantiateModelFactory() {
  static oops::ModelMaker<MODEL, oops::ModelInterface<MODEL, Model> > makerModelMpas_("MPAS");

  static oops::LinearModelMaker<MODEL, oops::LinearModelInterface<MODEL, Tlm> >
                                                                        makerTLMPAS_("MPASTLM");
}

}  // namespace soca

