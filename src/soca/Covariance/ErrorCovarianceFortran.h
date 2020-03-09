/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_COVARIANCE_ERRORCOVARIANCEFORTRAN_H_
#define SOCA_COVARIANCE_ERRORCOVARIANCEFORTRAN_H_

#include "soca/Fortran.h"
#include "soca/Geometry/Geometry.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_b_setup_f90(F90bmat &, const eckit::Configuration * const *,
                          const Geometry::Ftn * const &, const F90flds &,
                          const oops::Variables &);
    void soca_b_delete_f90(F90bmat &);
    void soca_b_mult_f90(const F90bmat &, const F90flds &,
                         const F90flds &);
    void soca_b_invmult_f90(const F90bmat &, const F90flds &, const F90flds &);
    void soca_b_randomize_f90(const F90bmat &, const F90flds &);
  }
}  // namespace soca
#endif  // SOCA_COVARIANCE_ERRORCOVARIANCEFORTRAN_H_
