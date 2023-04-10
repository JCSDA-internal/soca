/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "soca/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_vertconv_setup_f90(F90balopmat &,
                                 const eckit::Configuration * const *,
                                 const F90flds &,
                                 const F90geom &);
    void soca_vertconv_delete_f90(F90balopmat &);
    void soca_vertconv_mult_f90(const F90balopmat &, F90balopmat &,
                                const F90balopmat &);
    void soca_vertconv_multad_f90(const F90balopmat &, F90balopmat &,
                                  const F90balopmat &);
  }
}  // namespace soca
