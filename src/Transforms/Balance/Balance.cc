/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Transforms/Balance/Balance.h"

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
  Balance::Balance(const State & bkg,
         	 const State & traj,
		 const Geometry & geom,
	         const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_balance_setup_f90(keyFtnConfig_, &configc, traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  Balance::~Balance() {
    soca_balance_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Balance::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_balance_mult_f90(keyFtnConfig_,
    			 dxa.fields().toFortran(),
    			 dxm.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    // dxa = K^-1 dxm
    soca_balance_multinv_f90(keyFtnConfig_,
    			 dxm.fields().toFortran(),
    			 dxa.fields().toFortran());    
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    std::cout<<"keyFtnConfig_:"<<keyFtnConfig_<<std::endl;
    //dxa = dxm;
    soca_balance_multad_f90(keyFtnConfig_,
    			 dxm.fields().toFortran(),
    			 dxa.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void Balance::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    // dxm = (K^-1)^T dxa
    soca_balance_multinvad_f90(keyFtnConfig_,
    			 dxa.fields().toFortran(),
    			 dxm.fields().toFortran());
    
  }
  // -----------------------------------------------------------------------------
  void Balance::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
