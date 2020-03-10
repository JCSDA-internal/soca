/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_VERTCONV_VERTCONVFORTRAN_H_
#define SOCA_TRANSFORMS_VERTCONV_VERTCONVFORTRAN_H_

#include "soca/Fields/Fields.h"
#include "soca/Transforms/VertConv/VertConv.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_vertconv_setup_f90(VertConv::Ftn * &,
                                 const eckit::Configuration * const &,
                                 const Fields::Ftn * const &,
                                 const Fields::Ftn * const &);
    void soca_vertconv_delete_f90(VertConv::Ftn * const &);
    void soca_vertconv_mult_f90(const VertConv::Ftn * const &,
                                const Fields::Ftn * const &,
                                Fields::Ftn * const &);
    void soca_vertconv_multad_f90(const VertConv::Ftn * const &,
                                  const Fields::Ftn * const &,
                                  Fields::Ftn * const &);
  }
}  // namespace soca
#endif  // SOCA_TRANSFORMS_VERTCONV_VERTCONVFORTRAN_H_
