/*
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Fortran.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/Localization/Localization.h"

#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace soca {
// -----------------------------------------------------------------------------
Localization::Localization(const Geometry & resol,
                                       const eckit::Configuration & config) {
  const eckit::Configuration * configc = &config;
  soca_localization_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
}
// -----------------------------------------------------------------------------
Localization::~Localization() {
  soca_localization_delete_f90(keyFtnConfig_);
}
// -----------------------------------------------------------------------------
void Localization::multiply(Increment & dx) const {
  soca_localization_mult_f90(keyFtnConfig_, dx.fields().toFortran());
}
// -----------------------------------------------------------------------------
void Localization::print(std::ostream & os) const {
  os << "Localization::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace soca
