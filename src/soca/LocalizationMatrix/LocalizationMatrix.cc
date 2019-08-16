/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "soca/LocalizationMatrix/LocalizationMatrix.h"
#include "soca/Fortran.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace soca {
// -----------------------------------------------------------------------------
LocalizationMatrix::LocalizationMatrix(const Geometry & resol,
                                       const eckit::Configuration & config) {
  const eckit::Configuration * configc = &config;
  soca_localization_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
}
// -----------------------------------------------------------------------------
LocalizationMatrix::~LocalizationMatrix() {
  soca_localization_delete_f90(keyFtnConfig_);
}
// -----------------------------------------------------------------------------
void LocalizationMatrix::multiply(Increment & dx) const {
  soca_localization_mult_f90(keyFtnConfig_, dx.fields().toFortran());
}
// -----------------------------------------------------------------------------
void LocalizationMatrix::print(std::ostream & os) const {
  os << "LocalizationMatrix::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace soca
