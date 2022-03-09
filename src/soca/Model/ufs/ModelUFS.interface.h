/*
 * (C) Copyright 2020 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "soca/Fortran.h"
#include "soca/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace soca {

extern "C" {

  void soca_ufs_create_f90(F90model &, const eckit::Configuration &, const F90geom &);
  void soca_ufs_delete_f90(F90model &);


  void soca_ufs_initialize_f90(const F90model &, const F90state &, util::DateTime * const *);

  void soca_ufs_step_f90(const F90model &, const F90state &, util::DateTime * const *,
                            util::DateTime * const *);

  void soca_ufs_finalize_f90(const F90model &, const F90inc &, util::DateTime * const *);

}
// -----------------------------------------------------------------------------
}  // namespace soca
