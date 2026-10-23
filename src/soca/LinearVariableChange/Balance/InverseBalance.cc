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
#include "oops/util/Timer.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/LinearVariableChange/Balance/InverseBalance.h"
#include "soca/LinearVariableChange/Balance/BalanceFortran.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------

  static LinearVariableChangeMaker<InverseBalance>
                         makerLinearVariableChangeBalance_("Inverse BalanceSOCA");

  // -----------------------------------------------------------------------------
  InverseBalance::InverseBalance(const State & bkg,
                                 const State & traj,
                                 const Geometry & geom,
                                 const eckit::Configuration & conf) {
    oops::Log::trace() << "soca::InverseBalance::setup " << std::endl;
    util::Timer timer("soca::InverseBalance", "InverseBalance");

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
  InverseBalance::~InverseBalance() {
    oops::Log::trace() << "soca::InverseBalance::delete " << std::endl;
    soca_balance_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void InverseBalance::multiply(const Increment & dxm, Increment & dxa) const {
    // dxa = K^-1 dxm
    oops::Log::trace() << "soca::InverseBalance::multiply " << std::endl;
    util::Timer timer("soca::InverseBalance", "multiply");
    soca_balance_multinv_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void InverseBalance::multiplyAD(const Increment & dxa,
                           Increment & dxm) const {
    // dxm = (K^-1)^T dxa
    oops::Log::trace() << "soca::InverseBalance::multiplyAD " << std::endl;
    util::Timer timer("soca::InverseBalance", "multiplyAD");
    soca_balance_multinvad_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void InverseBalance::print(std::ostream & os) const {
    os << "SOCA linear change variable: inverse Balance";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
