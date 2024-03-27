/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FieldSetHelpers.h"

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
  const atlas::FunctionSpace & fs = geom_.functionSpace();
  auto v_norm = atlas::array::make_view<double, 2>(normHz_);

  // horizontal normalization
  for (atlas::Field field : fset) {
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t i = 0; i < fs.size(); i++) {
      for (size_t lvl = 0; lvl < field.shape(1); lvl++){
        view(i, lvl) *= v_norm(i,lvl);
      }
    }
  }

  // horizontal and vertical diffusion
  diffusion_->multiply(fset);

  // horizontal normalization
  for (atlas::Field field : fset) {
    auto view = atlas::array::make_view<double, 2>(field);
    for (size_t i = 0; i < fs.size(); i++) {
      for (size_t lvl = 0; lvl < field.shape(1); lvl++){
        view(i, lvl) *= v_norm(i,lvl);
      }
    }
  }
}

// --------------------------------------------------------------------------------------

void SaberDiffusion::read() {
  // eckit::LocalConfiguration conf = *params_.read.value();


  // ------------------------------------------------------------------------------------
  // TODO, actually read in information
  // ------------------------------------------------------------------------------------

  // generate some dummy scales... that should be read in
  atlas::Field hzScales = geom_.functionSpace().createField<double>(
    atlas::option::levels(25));
  auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);
  v_hzScales.assign(800e3);
  
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

  const int LEVELS = 25;
  // ------------------------------------------------------------------------------------
  // TODO read/generate scales
  // ------------------------------------------------------------------------------------
  read();

  // TODO save the diffusion parameters to be read directly next time

  // ------------------------------------------------------------------------------------
  // horizontal calibration
  // ------------------------------------------------------------------------------------
  const int randomizationIterations = 50000;
  const atlas::FunctionSpace & fs = geom_.functionSpace();
  
  // fields that are needed to keep a running variance calculation
  atlas::Field s = fs.createField<double>(atlas::option::levels(LEVELS));
  atlas::Field m = fs.createField<double>(atlas::option::levels(LEVELS));  
  auto v_s = atlas::array::make_view<double, 2>(s);
  auto v_m = atlas::array::make_view<double, 2>(m);
  v_s.assign(0.0);
  v_m.assign(0.0);

  // create a config for a random variable
  oops::Variables vars;  
  const std::string RANDVAR = "randVar";
  {
    eckit::LocalConfiguration varConf;
    eckit::LocalConfiguration varConf2;
    varConf2.set("levels", LEVELS);
    varConf.set(RANDVAR, varConf2);
    vars = oops::Variables(varConf, std::vector<std::string>{RANDVAR});
  }
  oops::FieldSet3D rand(util::DateTime(), geom_.comm());    

  // Perform multiple iterations of calculating the variance of the diffusion
  // operator when random vectors are supplied
  for (size_t itr = 1; itr <= randomizationIterations; itr++ ){
    // generate random vector, and apply sqrt of horizontal diffusion
    rand.randomInit(fs, vars);    
    auto v_rand = atlas::array::make_view<double, 2>(rand[RANDVAR]);
    diffusion_->multiplySqrt(rand);

    // keep track of the stats needed for a running variance calculation
    // (Welford 1962 algorithm)
    double new_m;
    for (size_t i = 0; i < fs.size(); i++) {
      for (size_t lvl = 0; lvl < LEVELS; lvl++){
        double m = v_m(i, lvl);
        double f = v_rand(i, lvl);
        new_m = m + (f - m) / itr;
        v_s(i, lvl) += (f - m)*(f - new_m);
        v_m(i, lvl) = new_m;
      }      
    }
  }

  // calculate final normalization constants
  normHz_ = fs.createField<double>(atlas::option::levels(LEVELS));
  auto v_norm = atlas::array::make_view<double, 2>(normHz_);
  v_norm.assign(1.0);
  for (size_t i = 0; i < fs.size(); i++) {
    for (size_t lvl = 0; lvl < LEVELS; lvl++){
      if(v_s(i, lvl) > 0.0) {
         v_norm(i, lvl) = 1.0 / sqrt( v_s(i, lvl) / (randomizationIterations-1));
      }
    }
  }
    

  // ------------------------------------------------------------------------------------
  // vertical calibration
  // ------------------------------------------------------------------------------------

}

// --------------------------------------------------------------------------------------

} // namespace soca
