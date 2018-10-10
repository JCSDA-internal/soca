/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Transforms/Kst/Kst.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "src/Increment/Increment.h"
#include "src/State/State.h"
#include "src/Fortran.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  Kst::Kst(const State & bkg,
	   const State & traj,
	   const Geometry & geom,	   
	   const eckit::Configuration & conf): traj_(traj) {
    const eckit::Configuration * configc = &conf;
    oops::Log::trace() << "soca::Kst::setup " << std::endl;
    soca_kst_setup_f90(keyFtnConfig_, &configc, traj_.fields().toFortran());
  }
  // -----------------------------------------------------------------------------
  Kst::~Kst() {
    oops::Log::trace() << "soca::Kst::delete " << std::endl;    
    soca_kst_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Kst::multiply(const Increment & dxa, Increment & dxm) const {
    // dxm = K dxa
    oops::Log::trace() << "soca::Kst::multiply " << std::endl;    
    soca_kst_mult_f90(dxa.fields().toFortran(),
		      dxm.fields().toFortran(),
		      traj_.fields().toFortran(),
		      keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Kst::multiplyInverse(const Increment & dxm, Increment & dxa) const {
    oops::Log::trace() << "soca::Kst::multiplyInverse " << std::endl;    
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void Kst::multiplyAD(const Increment & dxm, Increment & dxa) const {
    // dxa = K^T dxm
    oops::Log::trace() << "soca::Kst::multiplyAD " << std::endl;
    soca_kst_multad_f90(dxm.fields().toFortran(),
			dxa.fields().toFortran(),
			traj_.fields().toFortran(),
		        keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void Kst::multiplyInverseAD(const Increment & dxa, Increment & dxm) const {
    oops::Log::trace() << "soca::Kst::multiplyInverseAD " << std::endl;    
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void Kst::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
