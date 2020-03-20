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
#include "soca/State/State.h"
#include "soca/Transforms/Balance/Balance.h"
#include "soca/Transforms/Balance/BalanceFortran.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  Balance::Balance(const State & bkg,
                   const State & traj,
                   const Geometry & geom,
                   const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_balance_setup_f90(keyFtnConfig_, &configc, traj_.toFortran());
  }
  // -----------------------------------------------------------------------------
  Balance::~Balance() {
    soca_balance_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Balance::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_balance_mult_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    // dxa = K^-1 dxm
    soca_balance_multinv_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    // dxa = dxm;
    soca_balance_multad_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyInverseAD(const Increment & dxa,
                                  Increment & dxm) const {
    // dxm = (K^-1)^T dxa
    soca_balance_multinvad_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
