/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "soca/Fields/Fields.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Transforms/BkgErrGodas/BkgErrGodas.h"
#include "soca/Transforms/BkgErrGodas/BkgErrGodasFortran.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  BkgErrGodas::BkgErrGodas(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_bkgerrgodas_setup_f90(keyFtnConfig_,
                               &configc,
                               traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  BkgErrGodas::~BkgErrGodas() {
    soca_bkgerrgodas_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_bkgerrgodas_mult_f90(keyFtnConfig_,
                         dxa.fields().toFortran(),
                         dxm.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::multiplyInverse(const Increment & dxm,
                                    Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    soca_bkgerrgodas_mult_f90(keyFtnConfig_,
                         dxm.fields().toFortran(),
                         dxa.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::multiplyInverseAD(const Increment & dxa,
                                      Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
