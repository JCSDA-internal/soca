/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Traits.h"
#include "soca/LinearVariableChange/HorizFilt/HorizFilt.h"
#include "soca/LinearVariableChange/HorizFilt/HorizFiltFortran.h"


namespace soca {

  // -----------------------------------------------------------------------------

  static LinearVariableChangeMaker<HorizFilt>
               makerLinearVariableChangeHorizFilt_("HorizFiltSOCA");

  // -----------------------------------------------------------------------------
  HorizFilt::HorizFilt(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf):
    geom_(new Geometry(geom)),
    vars_(conf, "filter variables") {
    const eckit::Configuration * configc = &conf;

    // Interpolate trajectory to the geom resolution
    State traj_at_geomres(geom, traj);

    // Compute averaging weights
    soca_horizfilt_setup_f90(keyFtnConfig_,
                             &configc,
                             geom_->toFortran(),
                             traj_at_geomres.toFortran(),
                             vars_);

    // Get number of iterations
    niter_ = configc->getInt("niter");
  }
  // -----------------------------------------------------------------------------
  HorizFilt::~HorizFilt() {
    soca_horizfilt_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiply(const Increment & dxin, Increment & dxout) const {
    dxout = dxin;
    Increment dx_tmp(dxin);
    for (unsigned int iter = 0; iter < niter_; ++iter) {
      dx_tmp = dxout;
      soca_horizfilt_mult_f90(keyFtnConfig_,
                              dx_tmp.toFortran(),
                              dxout.toFortran(),
                              geom_->toFortran());
    }
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiplyInverse(const Increment & dxin, Increment & dxout)
    const {
    dxout = dxin;
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiplyAD(const Increment & dxin, Increment & dxout) const {
    dxout = dxin;
    Increment dx_tmp(dxin);
    for (unsigned int iter = 0; iter < niter_; ++iter) {
      dx_tmp = dxout;
      soca_horizfilt_multad_f90(keyFtnConfig_,
                                dx_tmp.toFortran(),
                                dxout.toFortran(),
                                geom_->toFortran());
    }
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiplyInverseAD(const Increment & dxin, Increment & dxout)
    const {
    dxout = dxin;
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::print(std::ostream & os) const {
    os << "SOCA linear change variable: HorizFilt";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
