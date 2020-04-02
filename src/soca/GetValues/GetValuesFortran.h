/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GETVALUES_GETVALUESFORTRAN_H_
#define SOCA_GETVALUES_GETVALUESFORTRAN_H_

#include "soca/Fortran.h"

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
  void soca_getvalues_create_f90(F90getvalues &, const F90geom &, const F90locs &);
  void soca_getvalues_delete_f90(F90getvalues &);
  void soca_getvalues_fill_geovals_f90(const F90getvalues &, const F90geom &, const F90state &,
                                       const util::DateTime **, const util::DateTime **,
                                       const F90locs &, const F90goms &);

};  // extern "C"

// -------------------------------------------------------------------------------------------------

}  // namespace soca
#endif // SOCA_GETVALUES_GETVALUESFORTRAN_H_
