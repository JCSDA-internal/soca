/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "soca/SaberBlocks/Diffusion/SaberDiffusion.h"
#include "soca/Utils/Diffusion.h"

namespace soca {

// --------------------------------------------------------------------------------------

static saber::SaberCentralBlockMaker<SaberDiffusion> makerSaberDiffusion_("DIFFUSION");

// --------------------------------------------------------------------------------------

SaberDiffusion::SaberDiffusion(
    const oops::GeometryData & geometryData,
    const oops::Variables & centralVars,
    const eckit::Configuration & covarConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberCentralBlockBase(params, xb.validTime()),
    geom_(geometryData),
    params_(params),
    diffusion_(new Diffusion(geometryData))
{

}

// --------------------------------------------------------------------------------------

void SaberDiffusion::randomize(oops::FieldSet3D & fset) const {
  throw eckit::NotImplemented("randomize not implemented yet for the block " + 
    this->blockName(),Here());
}

// --------------------------------------------------------------------------------------

void SaberDiffusion::multiply(oops::FieldSet3D & fset) const {
  diffusion_->multiply(fset);  
}

// --------------------------------------------------------------------------------------

void SaberDiffusion::read() {
  eckit::LocalConfiguration conf = *params_.read.value();


  // ------------------------------------------------------------------------------------
  // TODO, actually read in information
  // ------------------------------------------------------------------------------------

  // generate some dummy scales... that should be read in
  atlas::Field hzScales = geom_.functionSpace().createField<double>(
    atlas::option::levels(25));
  auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);
  v_hzScales.assign(500e3);
  
  // create vt scales
  atlas::Field vtScales = geom_.functionSpace().createField<double>(
    atlas::option::levels(25));
  auto v_vtScales = atlas::array::make_view<double, 2>(vtScales);
  v_vtScales.assign(5.0);

  // apply masking
  auto v_mask = atlas::array::make_view<double, 2>(geom_.getField("interp_mask"));
  for (size_t i = 0; i < hzScales.shape(0); i++) {
    if (!v_mask(i,0)) {
      for(size_t level = 0; level < hzScales.shape(1); level++) {
        v_hzScales(i, level) = 0.0;
        v_vtScales(i, level) = 0.0;
      }
    }
  }

  diffusion_->setScales(hzScales, vtScales);
}

// --------------------------------------------------------------------------------------

void SaberDiffusion::directCalibration(const oops::FieldSets &) {
  eckit::LocalConfiguration conf = *params_.calibration.value();

  // ------------------------------------------------------------------------------------
  // TODO read/generate scales
  // ------------------------------------------------------------------------------------
  // generate some dummy scales... that should be read in
  atlas::Field hzScales = geom_.functionSpace().createField<double>(
    atlas::option::levels(25));
  auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);
  v_hzScales.assign(500e3);
  
  // create vt scales
  atlas::Field vtScales = geom_.functionSpace().createField<double>(
    atlas::option::levels(25));
  auto v_vtScales = atlas::array::make_view<double, 2>(vtScales);
  v_vtScales.assign(5.0);

  // apply masking
  auto v_mask = atlas::array::make_view<double, 2>(geom_.getField("interp_mask"));
  for (size_t i = 0; i < hzScales.shape(0); i++) {
    if (!v_mask(i,0)) {
      for(size_t level = 0; level < hzScales.shape(1); level++) {
        v_hzScales(i, level) = 0.0;
        v_vtScales(i, level) = 0.0;
      }
    }
  }

  diffusion_->setScales(hzScales, vtScales);

  // TODO save the diffusion parameters to be read directly next time

  // ------------------------------------------------------------------------------------
  // horizontal calibration
  // ------------------------------------------------------------------------------------
  

  // ------------------------------------------------------------------------------------
  // vertical calibration
  // ------------------------------------------------------------------------------------

}

// --------------------------------------------------------------------------------------

} // namespace soca
