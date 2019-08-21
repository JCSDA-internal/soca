/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Transforms/HorizFilt/HorizFilt.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "eckit/config/Configuration.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Transforms/HorizFilt/HorizFiltFortran.h"

namespace soca {
  // -----------------------------------------------------------------------------
  HorizFilt::HorizFilt(const State & bkg,
                 const State & traj,
                 const Geometry & geom,
                 const eckit::Configuration & conf):
    geom_(new Geometry(geom)), vars_(conf) {
    const eckit::Configuration * configc = &conf;
    const eckit::Configuration * confvars = &vars_.toFortran();
    soca_horizfilt_setup_f90(keyFtnConfig_,
			     &configc,
			     geom_->toFortran(),
			     &confvars);
  }
  // -----------------------------------------------------------------------------
  HorizFilt::~HorizFilt() {
    soca_horizfilt_delete_f90(keyFtnConfig_);
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiply(const Increment & dxa, Increment & dxm) const {
    //dxm = dxa;
    soca_horizfilt_mult_f90(keyFtnConfig_,
    			    dxa.fields().toFortran(),
    			    dxm.fields().toFortran(),
    			    geom_->toFortran());

  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiplyInverse(const Increment & dxm, Increment & dxa)
    const {
    dxa = dxm;
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiplyAD(const Increment & dxm, Increment & dxa) const {
    //dxa = dxm;
    soca_horizfilt_multad_f90(keyFtnConfig_,
    			      dxm.fields().toFortran(),
   			      dxa.fields().toFortran(),
    			      geom_->toFortran());
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::multiplyInverseAD(const Increment & dxa, Increment & dxm)
    const {
    dxm = dxa;
  }
  // -----------------------------------------------------------------------------
  void HorizFilt::print(std::ostream & os) const {
    os << "SOCA change variable";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
