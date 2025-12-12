/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#pragma once

#include "soca/Fortran.h"

namespace soca {
  extern "C" {
    void soca_soca2cice_setup_f90(F90varchange &,
                                  const eckit::Configuration &,
                                  const F90geom &);
    void soca_soca2cice_changevar_f90(const F90varchange &,
                                      const F90geom &,
                                      const F90flds &,
                                      F90flds &);
  }
}  // namespace soca
