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
#include "soca/Transforms/BkgErrFilt/BkgErrFilt.h"
#include "soca/Transforms/BkgErrFilt/BkgErrFiltFortran.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  BkgErrFilt::BkgErrFilt(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_bkgerrfilt_setup_f90(ftn_,
                              &configc,
                              traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  BkgErrFilt::~BkgErrFilt() {
    soca_bkgerrfilt_delete_f90(ftn_);
  }
  // -----------------------------------------------------------------------------
  void BkgErrFilt::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_bkgerrfilt_mult_f90(ftn_,
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
    soca_bkgerrfilt_mult_f90(ftn_,
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
