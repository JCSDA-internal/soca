/*
 * (C) Copyright 2021-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Traits.h"
#include "soca/LinearVariableChange/LinearModel2GeoVaLs/LinearModel2GeoVaLs.h"
#include "soca/LinearVariableChange/LinearModel2GeoVaLs/LinearModel2GeoVaLsFortran.h"

namespace soca {

static LinearVariableChangeMaker<LinearModel2GeoVaLs>
       makerLinearVariableChangeModel2GeoVaLs_("LinearModel2GeoVaLs");
static LinearVariableChangeMaker<LinearModel2GeoVaLs>
       makerLinearVariableChangeModel2GeoVaLsDefault_("default");

// -----------------------------------------------------------------------------

LinearModel2GeoVaLs::LinearModel2GeoVaLs(const State & bg, const State &fg,
                                        const Geometry &geom,
                                        const eckit::Configuration &conf)
  : geom_(new Geometry(geom)) {
}

// -----------------------------------------------------------------------------

LinearModel2GeoVaLs::~LinearModel2GeoVaLs() {
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiply(const Increment &dxin,
                                         Increment &dxout) const {
  soca_model2geovals_linear_changevar_f90(geom_->toFortran(),
                                          dxin.toFortran(), dxout.toFortran());
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiplyInverse(const Increment &dxin,
                                                Increment &dxout) const {
  multiply(dxin, dxout);
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiplyAD(const Increment &dxin,
                                           Increment &dxout) const {
  soca_model2geovals_linear_changevarAD_f90(geom_->toFortran(),
                                            dxin.toFortran(),
                                            dxout.toFortran());
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiplyInverseAD(const Increment &dxin,
                                                  Increment &dxout) const {
  multiplyAD(dxin, dxout);
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::print(std::ostream & os) const {
  os << "SOCA linear change variable: LinearModel2GeoVaLs";
}

// -----------------------------------------------------------------------------

}  // namespace soca
