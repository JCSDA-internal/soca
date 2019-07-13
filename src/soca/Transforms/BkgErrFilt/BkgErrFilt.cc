/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Transforms/BkgErrFilt/BkgErrFilt.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "BkgErrFiltFortran.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Geometry/Geometry.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  BkgErrFilt::BkgErrFilt(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_bkgerrfilt_setup_f90(keyFtnConfig_,
                              &configc,
                              traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  BkgErrFilt::~BkgErrFilt() {
    soca_bkgerrfilt_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_bkgerrfilt_mult_f90(keyFtnConfig_,
                         dxa.fields().toFortran(),
                         dxm.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiplyInverse(const Increment & dxm,
                                   Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    soca_bkgerrfilt_mult_f90(keyFtnConfig_,
                         dxm.fields().toFortran(),
                         dxa.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiplyInverseAD(const Increment & dxa,
                                     Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
