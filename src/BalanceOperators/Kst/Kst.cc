/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/BalanceOperators/Kst/Kst.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "src/Increment/Increment.h"
#include "src/State/State.h"
#include "src/Fortran.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  Kst::Kst(const State & bkg,
	   const State & traj,
	   const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_kst_setup_f90(keyFtnConfig_, &configc);
  }
  // -----------------------------------------------------------------------------
  Kst::~Kst() {
    soca_kst_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Kst::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_kst_mult_f90(dxa.fields().toFortran(),
		      dxm.fields().toFortran(),
		      traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void Kst::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void Kst::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm  
    soca_kst_multad_f90(dxm.fields().toFortran(),
			dxa.fields().toFortran(),
			traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void Kst::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void Kst::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
