/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_MOM6SOLO_MODELFORTRAN_H_
#define SOCA_MODEL_MOM6SOLO_MODELFORTRAN_H_

#include "soca/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_model_setup_f90(const eckit::Configuration * const *,
                              const F90geom &, F90model &);
    void soca_model_delete_f90(F90model &);
    void soca_model_init_f90(const F90model &, const F90flds &);
    void soca_model_finalize_f90(const F90model &, const F90flds &);
    void soca_model_propagate_f90(const F90model &, const F90flds &,
                                  util::DateTime * const *);
  }
}  // namespace soca

#endif  // SOCA_MODEL_MOM6SOLO_MODELFORTRAN_H_
