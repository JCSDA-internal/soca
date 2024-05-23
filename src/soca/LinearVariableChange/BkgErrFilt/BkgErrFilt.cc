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
#include "soca/Traits.h"
#include "soca/LinearVariableChange/BkgErrFilt/BkgErrFilt.h"
#include "soca/LinearVariableChange/BkgErrFilt/BkgErrFiltFortran.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------

  static LinearVariableChangeMaker<BkgErrFilt>
                         makerLinearVariableChangeBkgErrFilt_("BkgErrFILT");

  // -----------------------------------------------------------------------------
  BkgErrFilt::BkgErrFilt(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf) {
    util::Timer timer("soca::BkgErrFilt", "BkgErrFilt");
    const eckit::Configuration * configc = &conf;

    // Interpolate background to the geom resolution
    State bkg_at_geomres(geom, bkg);

    // Setup background error filter
    soca_bkgerrfilt_setup_f90(keyFtnConfig_,
                              &configc,
                              bkg_at_geomres.toFortran(),
                              geom.toFortran());
  }
  // -----------------------------------------------------------------------------
  BkgErrFilt::~BkgErrFilt() {
    util::Timer timer("soca::BkgErrFilt", "~BkgErrFilt");
    soca_bkgerrfilt_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiply(const Increment & dxa, Increment & dxm) const {
    util::Timer timer("soca::BkgErrFilt", "multiply");
    // dxm = K dxa
    soca_bkgerrfilt_mult_f90(keyFtnConfig_, dxa.toFortran(), dxm.toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiplyInverse(const Increment & dxm,
                                   Increment & dxa) const {
    util::Timer timer("soca::BkgErrFilt", "multiplyInverse");
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiplyAD(const Increment & dxm, Increment & dxa) const {
    util::Timer timer("soca::BkgErrFilt", "multiplyAD");
    // dxa = K^T dxm
    soca_bkgerrfilt_mult_f90(keyFtnConfig_, dxm.toFortran(), dxa.toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiplyInverseAD(const Increment & dxa,
                                     Increment & dxm) const {
    util::Timer timer("soca::BkgErrFilt", "multiplyInverseAD");
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::print(std::ostream & os) const {
    os << "SOCA linear change variable: BkgErrFilt";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
