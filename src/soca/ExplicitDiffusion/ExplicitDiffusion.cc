/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/ExplicitDiffusion/ExplicitDiffusion.h"

namespace soca {

// --------------------------------------------------------------------------------------

static saber::SaberCentralBlockMaker<ExplicitDiffusion> makerExplicitDiffusion_("EXPLICIT_DIFFUSION");

// --------------------------------------------------------------------------------------

ExplicitDiffusion::ExplicitDiffusion(
    const oops::GeometryData & geometryData,
    const oops::Variables & centralVars,
    const eckit::Configuration & covarConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberCentralBlockBase(params)
{
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::randomize(atlas::FieldSet &) const {

}

// --------------------------------------------------------------------------------------
  
void ExplicitDiffusion::multiply(atlas::FieldSet &) const {

}

void ExplicitDiffusion::directCalibration(const std::vector<atlas::FieldSet> &) {
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::read() {
}

// --------------------------------------------------------------------------------------

void ExplicitDiffusion::print(std::ostream &) const {

}

// --------------------------------------------------------------------------------------

}  // namespace soca