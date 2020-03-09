/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_BKGERRFILT_BKGERRFILTFORTRAN_H_
#define SOCA_TRANSFORMS_BKGERRFILT_BKGERRFILTFORTRAN_H_

#include "soca/Fields/Fields.h"
#include "soca/Transforms/BkgErrFilt/BkgErrFilt.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_bkgerrfilt_setup_f90(BkgErrFilt::Ftn * &,
                                   const eckit::Configuration * const *,
                                   const Fields::Ftn * const &);
    void soca_bkgerrfilt_delete_f90(BkgErrFilt::Ftn * &);
    void soca_bkgerrfilt_mult_f90(const BkgErrFilt::Ftn * const &,
                                  const Fields::Ftn * const &,
                                  const Fields::Ftn * const &);
    void soca_bkgerrfilt_multad_f90(const BkgErrFilt::Ftn * const &,
                                    const Fields::Ftn * const &,
                                    const Fields::Ftn * const &);
  }
}  // namespace soca

#endif  // SOCA_TRANSFORMS_BKGERRFILT_BKGERRFILTFORTRAN_H_
