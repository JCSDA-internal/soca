/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_STATE_STATEFORTRAN_H_
#define SOCA_STATE_STATEFORTRAN_H_

#include "soca/Fortran.h"

#include "oops/base/Variables.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace soca {

  extern "C" {
    void soca_state_create_f90(F90flds &, const F90geom &,
                               const oops::Variables &);
    void soca_state_delete_f90(F90flds &);
    void soca_state_copy_f90(const F90flds &, const F90flds &);
    void soca_state_zero_f90(const F90flds &);
    void soca_state_axpy_f90(const F90flds &, const double &, const F90flds &);
    void soca_state_add_incr_f90(const F90flds &, const F90flds &);
    void soca_state_read_file_f90(const F90flds &,
                                  const eckit::Configuration * const &,
                                  util::DateTime * const *);
    void soca_state_write_file_f90(const F90flds &,
                                   const eckit::Configuration * const &,
                                   const util::DateTime * const *);
    void soca_state_rotate2grid_f90(const F90flds &,
                                    const oops::Variables &,
                                    const oops::Variables &);
    void soca_state_rotate2north_f90(const F90flds &,
                                     const oops::Variables &,
                                     const oops::Variables &);
    void soca_state_logtrans_f90(const F90flds &, const oops::Variables &);
    void soca_state_expontrans_f90(const F90flds &, const oops::Variables &);
    void soca_state_gpnorm_f90(const F90flds &, const int &, double &);
    void soca_state_sizes_f90(const F90flds &, int &,
                              int &, int &, int &);
    void soca_state_rms_f90(const F90flds &, double &);
    void soca_state_change_resol_f90(const F90flds &, const F90flds &);
    void soca_state_serial_size_f90(const F90flds &,
                                    const F90geom &,
                                    size_t &);
    void soca_state_serialize_f90(const F90flds &,
                                  const F90geom &,
                                  const size_t &,
                                  double[]);
    void soca_state_deserialize_f90(const F90flds &,
                                    const F90geom &,
                                    const size_t &,
                                    const double[],
                                    size_t &);
    void soca_state_analytic_f90(const F90flds &,
                                 const eckit::Configuration * const &,
                                 util::DateTime * const *);
    void soca_state_has_fields_f90(const F90flds &, const oops::Variables &, bool &);
    void soca_state_update_fields_f90(const F90flds &, const oops::Variables &);
  }
}  // namespace soca
#endif  // SOCA_STATE_STATEFORTRAN_H_
