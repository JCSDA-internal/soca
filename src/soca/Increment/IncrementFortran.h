/*
 * (C) Copyright 2020-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_INCREMENT_INCREMENTFORTRAN_H_
#define SOCA_INCREMENT_INCREMENTFORTRAN_H_

#include "oops/base/Variables.h"
#include "soca/Fortran.h"

// Forward declarations
namespace atlas {
  namespace field {
    class FieldSetImpl;
  }
}
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace soca {

  extern "C" {
    void soca_increment_create_f90(F90flds &, const F90geom &,
                               const oops::Variables &,
                               const atlas::field::FieldSetImpl *);
    void soca_increment_delete_f90(F90flds &);
    void soca_increment_copy_f90(const F90flds &, const F90flds &);
    void soca_increment_random_f90(const F90flds &);
    void soca_increment_dirac_f90(const F90flds &,
                              const eckit::Configuration * const &);
    void soca_increment_read_file_f90(const F90flds &,
                                  const eckit::Configuration * const &,
                                  util::DateTime * const *);
    void soca_increment_write_file_f90(const F90flds &,
                                   const eckit::Configuration * const &,
                                   const util::DateTime * const *);
    void soca_increment_update_fields_f90(F90flds &, const oops::Variables &);
    void soca_increment_horiz_scales_f90(F90flds &,
                                         const eckit::Configuration * const &);
    void soca_increment_vert_scales_f90(F90flds &, const double &);
  }
}  // namespace soca
#endif  // SOCA_INCREMENT_INCREMENTFORTRAN_H_
