/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef SOCA_VARIABLECHANGE_MODEL2GEOVALSFORTRAN_H_
#define SOCA_VARIABLECHANGE_MODEL2GEOVALSFORTRAN_H_

#include "soca/Fortran.h"

namespace soca {
  extern "C" {
    void soca_model2geovals_changevar_f90(const F90geom &, const F90flds &, F90flds &);
    void soca_model2geovals_linear_changevar_f90(const F90geom &, const F90flds &, F90flds &);
    void soca_model2geovals_linear_changevarAD_f90(const F90geom &, const F90flds &, F90flds &);
  }
}

#endif  // SOCA_VARIABLECHANGE_MODEL2GEOVALSFORTRAN_H_/
