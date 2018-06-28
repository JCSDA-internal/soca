/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/Covariance/ErrorCovariance.h"

#include <cmath>

#include "oops/util/Logger.h"
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
				   const eckit::Configuration & conf, const State & bkg) {
    //bkg: Background state, invariant wrt outer-loop. 
    time_ = util::DateTime(conf.getString("date"));
    const eckit::Configuration * configc = &conf;
    soca_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran(), bkg.fields().toFortran());
    Log::trace() << "ErrorCovariance created" << std::endl;
  }

  // -----------------------------------------------------------------------------

  ErrorCovariance::~ErrorCovariance() {
    soca_b_delete_f90(keyFtnConfig_);
    Log::trace() << "ErrorCovariance destructed" << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::linearize(const State & traj, const Geometry & resol) {
    geom_.reset(new Geometry(resol));
  //traj: Trajectory used for the linearization of the balance operators. Changes at each outer-loops.
    traj_.reset(new State(traj));
    soca_b_linearize_f90(traj.fields().toFortran(), resol.toFortran());
    Log::trace() << "Trajectory for ErrorCovariance" << std::endl;    
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::multiply(const Increment & dxin, Increment & dxout)
    const {
    Log::debug() << std::endl << "------- dxin" << dxin << std::endl;
    Log::debug() << std::endl <<"------ dxout" << dxout << std::endl;
    Log::debug() << std::endl <<"------ traj ---- :" << std::endl;
    soca_b_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
     		    dxout.fields().toFortran(), traj_->fields().toFortran());
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::inverseMultiply(const Increment & dxin, Increment & dxout) const {
    soca_b_invmult_f90(keyFtnConfig_, dxin.fields().toFortran(),		       
     			    dxout.fields().toFortran());
    std::cout << "inv mult" << std::endl;        
    //dxout = dxin;
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
