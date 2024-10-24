/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Traits.h"
#include "soca/LinearVariableChange/BkgErr/BkgErr.h"
#include "soca/LinearVariableChange/BkgErr/BkgErrFortran.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------

  static LinearVariableChangeMaker<BkgErr>
                         makerLinearVariableChangeBkgErr_("BkgErrSOCA");

  // -----------------------------------------------------------------------------
  BkgErr::BkgErr(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf) {
    util::Timer timer("soca::BkgErr", "BkgErr");
    const eckit::Configuration * configc = &conf;

    // Interpolate background to the geom resolution
    State bkg_at_geomres(geom, bkg);

    // Read/setup the diagonal of B
    soca_bkgerr_setup_f90(keyFtnConfig_,
                          &configc,
                          bkg_at_geomres.toFortran(),
                          geom.toFortran());
  }
  // -----------------------------------------------------------------------------
  BkgErr::~BkgErr() {
    util::Timer timer("soca::BkgErr", "~BkgErr");
    soca_bkgerr_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiply(const Increment & dxa, Increment & dxm) const {
    util::Timer timer("soca::BkgErr", "multiply");
    // dxm = K dxa
    soca_bkgerr_mult_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    util::Timer timer("soca::BkgErr", "multiply");
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiplyAD(const Increment & dxm, Increment & dxa) const {
    util::Timer timer("soca::BkgErr", "multiplyAD");
    // dxa = K^T dxm
    soca_bkgerr_mult_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErr::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    util::Timer timer("soca::BkgErr", "multiplyInverseAD");
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void BkgErr::print(std::ostream & os) const {
    os << "SOCA linear change variable: BkgErr";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
