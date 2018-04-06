/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/Covariance/ErrorCovariance.h"

#include <cmath>

#include "util/Logger.h"
#include "model/Fields/Fields.h"
#include "model/Fortran.h"
#include "model/Geometry/Geometry.h"
#include "model/Increment/Increment.h"
#include "model/State/State.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"


using oops::Log;

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------

  ErrorCovariance::ErrorCovariance(const Geometry & resol, const oops::Variables &,
				   const eckit::Configuration & conf, const State &) {
    time_ = util::DateTime(conf.getString("date"));
    const eckit::Configuration * configc = &conf;
    soca_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
    Log::trace() << "ErrorCovariance created" << std::endl;
  }

  // -----------------------------------------------------------------------------

  ErrorCovariance::~ErrorCovariance() {
    soca_b_delete_f90(keyFtnConfig_);
    Log::trace() << "ErrorCovariance destructed" << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::linearize(const State &, const Geometry & resol) {
    geom_.reset(new Geometry(resol));
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::multiply(const Increment & dxin, Increment & dxout) const {
    std::cout << "mult" << std::endl;    
    dxout = dxin;
    Log::debug() << std::endl << "------- dxin" << dxin << std::endl;
    Log::debug() << std::endl <<"------ dxout" << dxout << std::endl;        
    //soca_b_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
    //   			 dxout.fields().toFortran());
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::inverseMultiply(const Increment & dxin, Increment & dxout) const {
    //soca_b_invmult_f90(keyFtnConfig_, dxin.fields().toFortran(),		       
    // 			    dxout.fields().toFortran());
    std::cout << "inv mult" << std::endl;        
    dxout = dxin;
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::randomize(Increment & dx) const {
    soca_b_randomize_f90(keyFtnConfig_, dx.fields().toFortran());
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::print(std::ostream & os) const {
    os << "ErrorCovariance::print not implemented";
  }

  // -----------------------------------------------------------------------------

}  // namespace soca
