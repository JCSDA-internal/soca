/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_MODELFORTRAN_H_
#define SOCA_MODEL_MODELFORTRAN_H_

#include "soca/Fields/Fields.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Model/Model.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_setup_f90(Model::Ftn * &,
                        const eckit::Configuration * const &,
                        const Geometry::Ftn * const &);
    void soca_delete_f90(Model::Ftn * const &);
    void soca_initialize_integration_f90(Model::Ftn * const &,
                                         Fields::Ftn * const &);
    void soca_finalize_integration_f90(Model::Ftn * const &,
                                       Fields::Ftn * const &);
    void soca_propagate_f90(Model::Ftn * const &,
                            Fields::Ftn * const &,
                            const util::DateTime * const &);
  }
}  // namespace soca

#endif  // SOCA_MODEL_MODELFORTRAN_H_
