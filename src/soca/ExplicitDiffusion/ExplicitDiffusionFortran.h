/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "soca/Fortran.h"

namespace soca {

typedef int F90explicitdiffusion;

extern "C" {
  void soca_explicitdiffusion_setup_f90( F90explicitdiffusion &, const F90geom &);
  void soca_explicitdiffusion_calibrate_f90( F90explicitdiffusion &);
}

}  // namespace soca