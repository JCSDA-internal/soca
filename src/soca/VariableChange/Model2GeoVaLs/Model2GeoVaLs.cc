/*
 * (C) Copyright 2021-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/VariableChange/Model2GeoVaLs/Model2GeoVaLs.h"
#include "soca/VariableChange/Model2GeoVaLs/Model2GeoVaLsFortran.h"

#include "oops/util/abor1_cpp.h"

namespace soca {

// -----------------------------------------------------------------------------

static VariableChangeMaker<Model2GeoVaLs>
                             makerVariableChangeModel2GeoVaLs_("Model2GeoVaLs");

static VariableChangeMaker<Model2GeoVaLs>
                                         makerVariableChangeDefault_("default");

// -----------------------------------------------------------------------------

Model2GeoVaLs::Model2GeoVaLs(const Geometry & geom,
                             const eckit::Configuration & conf)
  : geom_(new Geometry(geom)) {
}

// -----------------------------------------------------------------------------

Model2GeoVaLs::~Model2GeoVaLs() {}

// -----------------------------------------------------------------------------

void Model2GeoVaLs::changeVar(const State & xin, State & xout) const {
  soca_model2geovals_changevar_f90(geom_->toFortran(),
                                   xin.toFortran(), xout.toFortran());
}

// -----------------------------------------------------------------------------

void Model2GeoVaLs::changeVarInverse(const State &, State &) const {
  util::abor1_cpp("Model2GeoVaLs::changeVarInverse not implemented");
}

// -----------------------------------------------------------------------------

}  // namespace soca
