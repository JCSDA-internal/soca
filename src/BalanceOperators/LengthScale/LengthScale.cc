/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/BalanceOperators/LengthScale/LengthScale.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "src/Increment/Increment.h"
#include "src/State/State.h"
#include "src/Geometry/Geometry.h"
#include "src/Fortran.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  LengthScale::LengthScale(const State & bkg,
         	 const State & traj,
		 const Geometry & geom,
	         const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_lengthscale_setup_f90(keyFtnConfig_, &configc);
  }
  // -----------------------------------------------------------------------------
  LengthScale::~LengthScale() {
    soca_lengthscale_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void LengthScale::multiply(const Increment & dxa, Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void LengthScale::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void LengthScale::multiplyAD(const Increment & dxm, Increment & dxa) const {
    dxa = dxm;    
  }
  // -----------------------------------------------------------------------------
  void LengthScale::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void LengthScale::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
