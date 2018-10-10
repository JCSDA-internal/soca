/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Transforms/Ktc/Ktc.h"

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
  Ktc::Ktc(const State & bkg,
	   const State & traj,
	   const Geometry & geom,
	   const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    soca_ktc_setup_f90(keyFtnConfig_, &configc);
  }
  // -----------------------------------------------------------------------------
  Ktc::~Ktc() {
    soca_ktc_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Ktc::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    soca_ktc_mult_f90(dxa.fields().toFortran(),
		      dxm.fields().toFortran(),
		      traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void Ktc::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void Ktc::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm  
    soca_ktc_multad_f90(dxm.fields().toFortran(),
			dxa.fields().toFortran(),
			traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  void Ktc::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void Ktc::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
