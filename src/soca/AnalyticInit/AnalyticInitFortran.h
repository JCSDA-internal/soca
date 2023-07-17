/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_ANALYTICINIT_ANALYTICINITFORTRAN_H_
#define SOCA_ANALYTICINIT_ANALYTICINITFORTRAN_H_

#include "soca/Fortran.h"

// Forward declarations
namespace ufo {
  class SampledLocations;
}

namespace soca {

  extern "C" {
    void soca_analytic_geovals_f90(F90goms &, const ufo::SampledLocations &);
  }
}  // namespace soca
#endif  // SOCA_ANALYTICINIT_ANALYTICINITFORTRAN_H_
