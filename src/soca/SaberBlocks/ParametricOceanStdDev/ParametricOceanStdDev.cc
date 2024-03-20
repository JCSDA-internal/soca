/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "soca/SaberBlocks/ParametricOceanStdDev/ParametricOceanStdDev.h"

#include "saber/blocks/SaberOuterBlockBase.h"

namespace soca {

static saber::SaberOuterBlockMaker<ParametricOceanStdDev> makerParametricOceanStdDev_(
  "ParametricOceanStdDev");

// ------------------------------------------------------------------------------------------------

ParametricOceanStdDev::ParametricOceanStdDev(
    const oops::GeometryData & outerGeometryData,
    const oops::Variables & outerVars,
    const eckit::Configuration & covarConfig,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    innerVars_(outerVars),
    bkgErr_(xb.validTime(),xb.commGeom())
{
  const int levels = xb["hocn"].levels();

  auto v_lonlat = atlas::array::make_view<double, 2>(innerGeometryData_.functionSpace().lonlat());
  auto v_mask = atlas::array::make_view<double, 2>(innerGeometryData_.getField("interp_mask"));
  auto v_hocn = atlas::array::make_view<double, 2>(xb["hocn"]);
  auto v_tocn = atlas::array::make_view<double, 2>(xb["tocn"]);
  auto v_depth = atlas::array::make_view<double, 2>(xb["layer_depth"]);


  //*************************************************************************************
  // calculate T background error
  //*************************************************************************************
  atlas::Field tocn_err = innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(levels) |
    atlas::option::name("tocn"));
  bkgErr_.add(tocn_err);
  {
    // parameters from config
    auto minVal = params.tocn.value().min.value();
    auto maxVal = params.tocn.value().max.value();
    auto efold = params.tocn.value().efold.value();

    // output field
    auto v_tocn_err = atlas::array::make_view<double, 2>(tocn_err);  
    v_tocn_err.assign(0.0);
    
    // loop over all points
    std::vector<double> dtdz(levels);
    for (size_t i = 0; i < tocn_err.shape(0); i++) {
      // skip land
      if (v_mask(i, 0) == 0.0) continue;
      
      // calculate dt/dz
      // TODO this breaks with the thin layers, need to impose a min thickness
      for (size_t z = 1; z < levels-1; z++) {
        dtdz[z] = (v_tocn(i, z+1) - v_tocn(i, z-1)) / 
                  (v_hocn(i,z) + 0.5 * (v_hocn(i, z+1) + v_hocn(i, z-1)));
      }
      dtdz[0] = dtdz[1];
      dtdz[levels-1] = dtdz[levels-2];

      // calculate value as function of dt/dz, efolding scale, and min/max
      for (size_t z = 0; z < levels; z++) {
        // step 2: calc value from dT/dz
        auto val = abs(params.tocn.value().dz * dtdz[z]);

        // step 2: calc a minimum from the efolding scale, and min SST value
        double sstVal = 0.1; //TODO read this in
        auto localMin = minVal + (sstVal - minVal) * exp(-v_depth(i,z) / efold );        
       
        // step 3: min/max
        val = std::clamp(val, localMin,  maxVal);

        // step 4: smoothing
        // TODO

        // done
        v_tocn_err(i, z) = val;
      }
    }
  }

  //*************************************************************************************
  // calculate unbalanced S background error
  //*************************************************************************************
  bkgErr_.add(
    innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(levels) |
    atlas::option::name("socn")));

  //*************************************************************************************
  // calculate unbalanced SSH background error
  //*************************************************************************************
  atlas::Field ssh_err = innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(1) |
    atlas::option::name("ssh"));
  bkgErr_.add(ssh_err);
  {
    // parameter from input config
    auto minVal = params.ssh.value().min.value();
    auto maxVal = params.ssh.value().max.value();
  
    // output field
    auto v_ssh_err = atlas::array::make_view<double, 2>(ssh_err);
    v_ssh_err.assign(0.0);

    // loop over all points
    for (size_t i = 0; i < ssh_err.shape(0); i++) {
      // skip land
      if (v_mask(i, 0) == 0.0) continue;

      auto absLat = std::abs(v_lonlat(i, 1));
      if (absLat >= params.ssh.value().phi_ex.value()) {
        // In the extratropics, set to max value
        v_ssh_err(i, 0) = maxVal;
      } else {
        // otherwise, taper to min value towards the equator        
        constexpr double pi = 3.14159265358979323846;
        v_ssh_err(i, 0) = minVal + 0.5 * (maxVal - minVal)*(1 - std::cos(pi * absLat / params.ssh.value().phi_ex.value()));
      }     
    }
  }

  // optionally write the output
  if (params.output.value() != boost::none) {
    bkgErr_.write(*params.output.value());
  }

}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiply(oops::FieldSet3D & fset) const {
  
  for (atlas::Field field : fset) {
    if (! bkgErr_.has(field.name())) continue;

    auto v_field = atlas::array::make_view<double, 2>(field);
    auto v_mult = atlas::array::make_view<double, 2>(bkgErr_[field.name()]);

    for (size_t i = 0; i < field.shape(0); i++) {
      for (size_t z = 0; z < field.levels(); z++) {
        // TEMP, change this to *=
        v_field(i, z) = v_mult(i, z);
      }
    }
  }
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiplyAD(oops::FieldSet3D & fset) const
{
  multiply(fset);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::leftInverseMultiply(oops::FieldSet3D &) const
{
  ASSERT(1 == 3);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::print(std::ostream & os) const {
  os << classname();
}

// ------------------------------------------------------------------------------------------------
}  // namespace soca
