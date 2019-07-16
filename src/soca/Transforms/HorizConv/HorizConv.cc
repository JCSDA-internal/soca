/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Transforms/HorizConv/HorizConv.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Geometry/Geometry.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  HorizConv::HorizConv(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
		 const eckit::Configuration & conf):
    vars_(conf), cov2d_(geom, vars_, conf, bkg, traj){
    std::cout << "=============== HORIZCONV ======================" << std::endl;
  }
  // -----------------------------------------------------------------------------
  HorizConv::~HorizConv() {
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiply(const Increment & dxa, Increment & dxm) const {
    cov2d_.multiply(dxa, dxm);
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiplyAD(const Increment & dxm, Increment & dxa) const {
    cov2d_.multiply(dxm, dxa);
  }
  // -----------------------------------------------------------------------------
  void HorizConv::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void HorizConv::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
