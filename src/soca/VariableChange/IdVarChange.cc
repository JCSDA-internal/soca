/*
 * (C) Copyright 2021-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"

#include "oops/base/VariableChangeBase.h"
#include "oops/generic/IdVariableChange.h"

namespace soca {
  static oops::GenericVariableChangeMaker
         <Traits, oops::IdVariableChange<Traits> >
         makerVariableChangeId_("default");
}
