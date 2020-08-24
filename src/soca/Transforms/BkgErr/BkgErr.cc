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
#include "soca/Transforms/BkgErr/BkgErr.h"
#include "soca/Transforms/BkgErr/BkgErrFortran.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  BkgErr::BkgErr(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf) :
      bkgerr_variance_(initBkgErr(geom, conf)) {
    oops::Log::trace() << "soca::BkgErr::setup " << std::endl;
    const eckit::Configuration * configc = &conf;

    // Setup diagonal standard deviation of the background error
    soca_bkgerr_setup_f90(keyFtnConfig_,
                          &configc,
                          bkgerr_variance_.toFortran(),
                          geom.toFortran());
  }
  // -----------------------------------------------------------------------------
  BkgErr::~BkgErr() {
    soca_bkgerr_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = D dxa
    soca_bkgerr_mult_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = D dxm
    soca_bkgerr_mult_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void BkgErr::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
