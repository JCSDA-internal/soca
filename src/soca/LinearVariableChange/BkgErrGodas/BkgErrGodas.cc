/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/LinearVariableChange/BkgErrGodas/BkgErrGodas.h"
#include "soca/LinearVariableChange/BkgErrGodas/BkgErrGodasFortran.h"
#include "soca/State/State.h"
#include "soca/Traits.h"

#include "eckit/config/Configuration.h"

#include "oops/interface/LinearVariableChange.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  static oops::LinearVariableChangeMaker<Traits,
              oops::LinearVariableChange<Traits, BkgErrGodas> >
              makerLinearVariableChangeBkgErrGodas_("BkgErrGODAS");

  // -----------------------------------------------------------------------------
  BkgErrGodas::BkgErrGodas(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf) {
    oops::Log::trace() << "soca::BkgErrGodas::setup " << std::endl;
    const eckit::Configuration * configc = &conf;

    // Interpolate trajectory to the geom resolution
    State traj_at_geomres(geom, traj);

    // Initialize the parametric background error variance
    soca_bkgerrgodas_setup_f90(keyFtnConfig_,
                               &configc,
                               traj_at_geomres.toFortran(),
                               geom.toFortran());
  }
  // -----------------------------------------------------------------------------
  BkgErrGodas::~BkgErrGodas() {
    oops::Log::trace() << "soca::BkgErrGodas::delete " << std::endl;
    soca_bkgerrgodas_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    oops::Log::trace() << "soca::BkgErrGodas::multiply " << std::endl;
    soca_bkgerrgodas_mult_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::multiplyInverse(const Increment & dxm,
                                    Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void BkgErrGodas::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    oops::Log::trace() << "soca::BkgErrGodas::multiplyAD " << std::endl;
    soca_bkgerrgodas_mult_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
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
