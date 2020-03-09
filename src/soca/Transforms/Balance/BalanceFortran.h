/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_BALANCE_BALANCEFORTRAN_H_
#define SOCA_TRANSFORMS_BALANCE_BALANCEFORTRAN_H_

#include "soca/Fields/Fields.h"
#include "soca/Transforms/Balance/Balance.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_balance_setup_f90(Balance::Ftn * &,
                                const eckit::Configuration * const *,
                                const Fields::Ftn * const &);
    void soca_balance_delete_f90(Balance::Ftn * &);
    void soca_balance_mult_f90(Balance::Ftn * const &,
                               const Fields::Ftn * const &,
                               Fields::Ftn * &);
    void soca_balance_multinv_f90(Balance::Ftn * const &,
                                  const Fields::Ftn * const &,
                                  Fields::Ftn * &);
    void soca_balance_multad_f90(Balance::Ftn * const &,
                                 const Fields::Ftn * const &,
                                 Fields::Ftn * &);
    void soca_balance_multinvad_f90(Balance::Ftn * const &,
                                    const Fields::Ftn * const &,
                                    Fields::Ftn * &);
  }
}  // namespace soca
#endif  // SOCA_TRANSFORMS_BALANCE_BALANCEFORTRAN_H_
