/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Traits.h"
#include "soca/LinearVariableChange/Balance/Balance.h"
#include "soca/LinearVariableChange/Balance/BalanceFortran.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------

  static LinearVariableChangeMaker<Balance>
                         makerLinearVariableChangeBalance_("BalanceSOCA");

  // -----------------------------------------------------------------------------
  Balance::Balance(const State & bkg,
                   const State & traj,
                   const Geometry & geom,
                   const eckit::Configuration & conf) {
    oops::Log::trace() << "soca::Balance::setup " << std::endl;
    const eckit::Configuration * configc = &conf;

    // Interpolate trajectory to the geom resolution
    State traj_at_geomres(geom, traj);

    // Compute Jacobians of the balance wrt traj
    soca_balance_setup_f90(keyFtnConfig_,
                           &configc,
                           traj_at_geomres.toFortran(),
                           geom.toFortran());
  }
  // -----------------------------------------------------------------------------
  Balance::~Balance() {
    oops::Log::trace() << "soca::Balance::delete " << std::endl;
    soca_balance_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Balance::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    oops::Log::trace() << "soca::Balance::multiply " << std::endl;
    soca_balance_mult_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    // dxa = K^-1 dxm
    oops::Log::trace() << "soca::Balance::multiplyInverse " << std::endl;
    soca_balance_multinv_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    oops::Log::trace() << "soca::Balance::multiplyAD " << std::endl;
    soca_balance_multad_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyInverseAD(const Increment & dxa,
                                  Increment & dxm) const {
    // dxm = (K^-1)^T dxa
    oops::Log::trace() << "soca::Balance::multiplyInverseAD " << std::endl;
    soca_balance_multinvad_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::print(std::ostream & os) const {
    os << "SOCA linear change variable: Balance";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
