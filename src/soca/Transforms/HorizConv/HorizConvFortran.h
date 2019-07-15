/*
 * (C) Copyright 2017-2019  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_TRANSFORMS_HORIZCONV_HORIZCONVFORTRAN_H_
#define SRC_TRANSFORMS_HORIZCONV_HORIZCONVFORTRAN_H_

#include "src/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_horizconv_setup_f90(F90balopmat &,
                               const eckit::Configuration * const *,
                               const F90flds &);
    void soca_horizconv_delete_f90(F90balopmat &);
    void soca_horizconv_mult_f90(const F90balopmat &,
                              const F90flds &,
                                  F90flds &);
    void soca_horizconv_multad_f90(const F90balopmat,
                                F90flds &,
                                const F90flds &);
  }
}  // namespace soca

#endif  // SRC_TRANSFORMS_HORIZCONV_HORIZCONVFORTRAN_H_
