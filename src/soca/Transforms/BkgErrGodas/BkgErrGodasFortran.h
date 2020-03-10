/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_BKGERRGODAS_BKGERRGODASFORTRAN_H_
#define SOCA_TRANSFORMS_BKGERRGODAS_BKGERRGODASFORTRAN_H_

#include "soca/Fields/Fields.h"
#include "soca/Transforms/BkgErrGodas/BkgErrGodas.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_bkgerrgodas_setup_f90(BkgErrGodas::Ftn * &,
                                    const eckit::Configuration * const &,
                                    const Fields::Ftn * const &);
    void soca_bkgerrgodas_delete_f90(BkgErrGodas::Ftn * const &);
    void soca_bkgerrgodas_mult_f90(const BkgErrGodas::Ftn * const &,
                                   const Fields::Ftn * const &,
                                   Fields::Ftn * const &);
  }
}  // namespace soca
#endif  // SOCA_TRANSFORMS_BKGERRGODAS_BKGERRGODASFORTRAN_H_
