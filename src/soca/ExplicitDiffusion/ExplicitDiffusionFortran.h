/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "soca/Fortran.h"

namespace eckit {
  class Configuration;
}

namespace soca {

typedef int F90explicitdiffusion;

extern "C" {
  void soca_explicitdiffusion_setup_f90(F90explicitdiffusion &, const F90geom &,
                                        const eckit::Configuration * const &);
  void soca_explicitdiffusion_calibrate_f90(const F90explicitdiffusion &,
                                            const eckit::Configuration * const &);
  void soca_explicitdiffusion_multiply_f90(const F90explicitdiffusion &,
                                           const F90flds &);
  void soca_explicitdiffusion_writeparams_f90(const F90explicitdiffusion &,
                                              const eckit::Configuration * const &);
  void soca_explicitdiffusion_readparams_f90(const F90explicitdiffusion &,
                                             const eckit::Configuration * const &);
}

}  // namespace soca
