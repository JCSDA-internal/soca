/*
 * (C) Copyright 2017-2019  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_MODEL_MODELFORTRAN_H_
#define SRC_MODEL_MODELFORTRAN_H_

#include "src/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_setup_f90(const eckit::Configuration * const *,
                        const F90geom &, F90model &);
    void soca_delete_f90(F90model &);
    void soca_initialize_integration_f90(const F90model &, const F90flds &);
    void soca_finalize_integration_f90(const F90model &, const F90flds &);
    void soca_propagate_f90(const F90model &, const F90flds &,
                            util::DateTime * const *);
  }
}  // namespace soca

#endif  // SRC_MODEL_MODELFORTRAN_H_
