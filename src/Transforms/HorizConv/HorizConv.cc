/*
 * (C) Copyright 2017-2019  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Transforms/HorizConv/HorizConv.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "src/Increment/Increment.h"
#include "src/State/State.h"
#include "src/Fortran.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  HorizConv::HorizConv(const State & bkg,
                     const State & traj,
                     const Geometry & geom,
                     const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;

    soca_horizconv_setup_f90(keyFtnConfig_,
                            &configc,
                            bkg.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  HorizConv::~HorizConv() {
    soca_horizconv_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_horizconv_mult_f90(dxa.fields().toFortran(),
                      dxm.fields().toFortran(),
                      keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    soca_horizconv_multad_f90(dxm.fields().toFortran(),
                        dxa.fields().toFortran(),
                        keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiplyInverseAD(const Increment & dxa,
                                   Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void HorizConv::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
