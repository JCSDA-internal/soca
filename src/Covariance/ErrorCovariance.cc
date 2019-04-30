/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Covariance/ErrorCovariance.h"

#include <cmath>

#include "oops/util/Logger.h"
#include "src/Fields/Fields.h"
#include "src/Fortran.h"
#include "src/Geometry/Geometry.h"
#include "src/Increment/Increment.h"
#include "src/State/State.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/assimilation/GMRESR.h"


using oops::Log;

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------

  ErrorCovariance::ErrorCovariance(const Geometry & resol,
                                   const oops::Variables &,
                                   const eckit::Configuration & conf,
                                   const State & bkg,
                                   const State & traj) {
    // bkg: Background state, invariant wrt outer-loop.
    time_ = util::DateTime(conf.getString("date"));
    const eckit::Configuration * configc = &conf;
    soca_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran(),
                     bkg.fields().toFortran());
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
    Log::trace() << "Trajectory for ErrorCovariance" << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::multiply(const Increment & dxin, Increment & dxout)
    const {
    dxout = dxin;
    soca_b_mult_f90(keyFtnConfig_, dxin.fields().toFortran(),
                    dxout.fields().toFortran());
    Log::trace() << "ErrorCovariance multiply" << std::endl;
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::inverseMultiply(const Increment & dxin,
                                        Increment & dxout) const {
    // oops::IdentityMatrix<Increment> Id;
    // dxout.zero();
    // GMRESR(dxout, dxin, *this, Id, 10, 1.0e-6);
    dxout = dxin;
    Log::trace() << "ErrorCovariance inversemultiply" << std::endl;
  }

  // -----------------------------------------------------------------------------

  //  void ErrorCovariance::doRandomize(Increment & dx) const {
  void ErrorCovariance::randomize(Increment & dx) const {
    soca_b_randomize_f90(keyFtnConfig_, dx.fields().toFortran());
  }

  // -----------------------------------------------------------------------------

  void ErrorCovariance::print(std::ostream & os) const {
    os << "ErrorCovariance::print not implemented";
  }

  // -----------------------------------------------------------------------------

}  // namespace soca
