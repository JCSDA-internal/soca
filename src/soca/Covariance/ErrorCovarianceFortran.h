/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_COVARIANCE_ERRORCOVARIANCEFORTRAN_H_
#define SOCA_COVARIANCE_ERRORCOVARIANCEFORTRAN_H_

#include "soca/Covariance/ErrorCovariance.h"
#include "soca/Fields/Fields.h"
#include "soca/Geometry/Geometry.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_b_setup_f90(ErrorCovariance::Ftn * &,
                          const eckit::Configuration * const &,
                          const Geometry::Ftn * const &,
                          const Fields::Ftn * const &,
                          const oops::Variables &);
    void soca_b_delete_f90(ErrorCovariance::Ftn * const &);
    void soca_b_mult_f90(const ErrorCovariance::Ftn * const &,
                         const Fields::Ftn * const &,
                         Fields::Ftn * const &);
    void soca_b_randomize_f90(const ErrorCovariance::Ftn * const &,
                              Fields::Ftn * const &);
  }
}  // namespace soca
#endif  // SOCA_COVARIANCE_ERRORCOVARIANCEFORTRAN_H_
