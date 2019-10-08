/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_HORIZFILT_HORIZFILTFORTRAN_H_
#define SOCA_TRANSFORMS_HORIZFILT_HORIZFILTFORTRAN_H_

#include "soca/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_horizfilt_setup_f90(F90balopmat &,
                                  const eckit::Configuration * const *,
                                  const F90geom &,
                                  const F90flds &,
                                  const eckit::Configuration * const *);
    void soca_horizfilt_delete_f90(F90balopmat &);
    void soca_horizfilt_mult_f90(const F90balopmat &,
                                 const F90flds &,
                                 F90flds &,
                                 const F90geom &);
    void soca_horizfilt_multad_f90(const F90balopmat &,
                                   const F90flds &,
                                   F90flds &,
                                   const F90geom &);
  }
}  // namespace soca

#endif  // SOCA_TRANSFORMS_HORIZFILT_HORIZFILTFORTRAN_H_
