/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Timer.h"

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
  : saber::SaberCentralBlockBase(params, xb.validTime()), params_(params)
{
  util::Timer timer("soca::ExplicitDiffusion", "ExplicitDiffusion");

  // setup geometry
  geom_.reset(new Geometry(params_.geometry.value(), geometryData.comm()));

  // setup the fortran code
  eckit::LocalConfiguration conf = params_.toConfiguration();
  soca_explicitdiffusion_setup_f90(keyFortran_, geom_->toFortran(), &conf);

  vars_ = params.activeVars.value().get_value_or(centralVars);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::randomize(oops::FieldSet3D & fset) const {
  util::Timer timer("soca::ExplicitDiffusion", "randomize");

  // Create random increments
  fset.randomInit(geom_->functionSpace(), fset.variables());
  Increment dx(*geom_, vars_, util::DateTime());
  dx.fromFieldSet(fset.fieldSet());

  // apply square root of diffusion
  soca_explicitdiffusion_multiply_f90(keyFortran_, dx.toFortran(), true);

  // copy back to fieldset
  dx.toFieldSet(fset.fieldSet());
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::multiply(oops::FieldSet3D & fset) const {
  util::Timer timer("soca::ExplicitDiffusion", "multiply");

  Increment dx(*geom_, vars_, util::DateTime());
  dx.fromFieldSet(fset.fieldSet());

  soca_explicitdiffusion_multiply_f90(keyFortran_, dx.toFortran());

  // TODO(Travis) not the cleanest, but I don't care since ExplicitDiffusion is
  // being rewritten soon anyway. The fortran code of ExplicitDiffusion is not
  // leaving alone variables it should leave alone, so we have to copy fields
  // back
  atlas::FieldSet fs2;
  dx.toFieldSet(fs2);
  for (auto & f : fs2) {
    auto view = atlas::array::make_view<double, 2>(f);
    auto otherView = atlas::array::make_view<double, 2>(fset.fieldSet().field(f.name()));
    for (atlas::idx_t jnode = 0; jnode < f.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < f.shape(1); ++jlevel) {
          otherView(jnode, jlevel) = view(jnode, jlevel);
      }
    }
  }
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::directCalibration(const oops::FieldSets &) {
  util::Timer timer("soca::ExplicitDiffusion", "directCalibration");

  eckit::LocalConfiguration conf = *params_.calibration.value();
  soca_explicitdiffusion_calibrate_f90(keyFortran_, &conf);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::read() {
  util::Timer timer("soca::ExplicitDiffusion", "read");

  eckit::LocalConfiguration conf = *params_.read.value();
  soca_explicitdiffusion_readparams_f90(keyFortran_, &conf);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::write() const {
  util::Timer timer("soca::ExplicitDiffusion", "write");

  eckit::LocalConfiguration conf = *params_.calibration.value();
  soca_explicitdiffusion_writeparams_f90(keyFortran_, &conf);
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::print(std::ostream &) const {
}

// --------------------------------------------------------------------------------------

}  // namespace soca
