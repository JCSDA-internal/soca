/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/ExplicitDiffusion/ExplicitDiffusion.h"
#include "soca/ExplicitDiffusion/ExplicitDiffusionFortran.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"

namespace soca {

// --------------------------------------------------------------------------------------

static saber::SaberCentralBlockMaker<ExplicitDiffusion>
  makerExplicitDiffusion_("EXPLICIT_DIFFUSION");

// --------------------------------------------------------------------------------------

ExplicitDiffusion::ExplicitDiffusion(
    const oops::GeometryData & geometryData,
    const oops::Variables & centralVars,
    const eckit::Configuration & covarConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberCentralBlockBase(params), params_(params)
{
  // setup geometry
  geom_.reset(new Geometry(params_.geometry.value(), geometryData.comm()));

  // setup the fortran code
  soca_explicitdiffusion_setup_f90(keyFortran_, geom_->toFortran());

  vars_ = params.activeVars.value().get_value_or(centralVars);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::randomize(atlas::FieldSet &) const {
  throw eckit::NotImplemented("read not implemented yet for ExplictDiffusion");
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::multiply(atlas::FieldSet & fset) const {
  Increment dx(*geom_, vars_, util::DateTime());
  dx.fromFieldSet(fset);

  soca_explicitdiffusion_multiply_f90(keyFortran_, dx.toFortran());

  dx.toFieldSet(fset);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::directCalibration(const std::vector<atlas::FieldSet> &) {
  eckit::LocalConfiguration conf = (*params_.calibration.value()).toConfiguration();
  soca_explicitdiffusion_calibrate_f90(keyFortran_, &conf);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::read() {
  eckit::LocalConfiguration conf = (*params_.read.value()).toConfiguration();
  soca_explicitdiffusion_readparams_f90(keyFortran_, &conf);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::write() const {
  eckit::LocalConfiguration conf =
    (*(*params_.calibration.value()).write.value()).toConfiguration();
  soca_explicitdiffusion_writeparams_f90(keyFortran_, &conf);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::print(std::ostream &) const {
}

// --------------------------------------------------------------------------------------

}  // namespace soca
