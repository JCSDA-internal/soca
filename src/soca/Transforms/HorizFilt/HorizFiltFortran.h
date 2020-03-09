/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_HORIZFILT_HORIZFILTFORTRAN_H_
#define SOCA_TRANSFORMS_HORIZFILT_HORIZFILTFORTRAN_H_

#include "soca/Fields/Fields.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Transforms/HorizFilt/HorizFilt.h"

#include "oops/base/Variables.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_horizfilt_setup_f90(HorizFilt::Ftn * &,
                                  const eckit::Configuration * const *,
                                  const Geometry::Ftn * const &,
                                  const Fields::Ftn * const &,
                                  const oops::Variables &);
    void soca_horizfilt_delete_f90(HorizFilt::Ftn * &);
    void soca_horizfilt_mult_f90(const HorizFilt::Ftn * const &,
                                 const Fields::Ftn * const &,
                                 const Fields::Ftn * const &);
    void soca_horizfilt_multad_f90(const HorizFilt::Ftn * const &,
                                   const Fields::Ftn * const &,
                                   const Fields::Ftn * const &);
  }
}  // namespace soca

#endif  // SOCA_TRANSFORMS_HORIZFILT_HORIZFILTFORTRAN_H_
