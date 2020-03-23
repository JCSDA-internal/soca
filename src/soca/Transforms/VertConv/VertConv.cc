/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Transforms/VertConv/VertConv.h"
#include "soca/Transforms/VertConv/VertConvFortran.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  VertConv::VertConv(const State & bkg,
                     const State & traj,
                     const Geometry & geom,
                     const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    oops::Log::trace() << "soca::VertConv::setup " << std::endl;
    soca_vertconv_setup_f90(keyFtnConfig_,
                            &configc,
                            traj_.toFortran(),
                            bkg.toFortran());
  }
  // -----------------------------------------------------------------------------
  VertConv::~VertConv() {
    oops::Log::trace() << "soca::VertConv::delete " << std::endl;
    soca_vertconv_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void VertConv::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    oops::Log::trace() << "soca::VertConv::multiply " << std::endl;
    soca_vertconv_mult_f90(dxa.toFortran(), dxm.toFortran(), keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void VertConv::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    oops::Log::trace() << "soca::VertConv::multiplyInverse " << std::endl;
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void VertConv::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    oops::Log::trace() << "soca::VertConv::multiplyAD " << std::endl;
    soca_vertconv_multad_f90(dxm.toFortran(), dxa.toFortran(), keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void VertConv::multiplyInverseAD(const Increment & dxa,
                                   Increment & dxm) const {
    oops::Log::trace() << "soca::VertConv::multiplyInverseAD " << std::endl;
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void VertConv::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
