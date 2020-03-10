/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_BKGERR_BKGERRFORTRAN_H_
#define SOCA_TRANSFORMS_BKGERR_BKGERRFORTRAN_H_

#include "soca/Fields/Fields.h"
#include "soca/Transforms/BkgErr/BkgErr.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

// -----------------------------------------------------------------------------

namespace soca {

  extern "C" {
    void soca_bkgerr_setup_f90(BkgErr::Ftn * &,
                               const eckit::Configuration * const &,
                               const Fields::Ftn * const &);
    void soca_bkgerr_delete_f90(BkgErr::Ftn * const &);
    void soca_bkgerr_mult_f90(const BkgErr::Ftn * const &,
                              const Fields::Ftn * const &,
                              Fields::Ftn * const &);
  }
}  // namespace soca

#endif  // SOCA_TRANSFORMS_BKGERR_BKGERRFORTRAN_H_
