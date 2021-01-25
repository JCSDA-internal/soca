/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_LOCALIZATION_LOCALIZATIONFORTRAN_H_
#define SOCA_LOCALIZATION_LOCALIZATIONFORTRAN_H_

#include "soca/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_localization_setup_f90(F90lclz &,
                                     const eckit::Configuration * const *,
                                     const F90geom &);
    void soca_localization_delete_f90(F90lclz &);
    void soca_localization_randomize_f90(const F90lclz &, const F90flds &);
    void soca_localization_mult_f90(const F90lclz &, const F90flds &);
  }
}  // namespace soca
#endif  // SOCA_LOCALIZATION_LOCALIZATIONFORTRAN_H_
