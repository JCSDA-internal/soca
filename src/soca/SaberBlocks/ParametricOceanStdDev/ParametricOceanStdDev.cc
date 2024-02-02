/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

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
    innerVars_(outerVars)
{
  const int levels = xb["hocn"].levels();

  auto v_mask = atlas::array::make_view<double, 2>(innerGeometryData_.getField("interp_mask"));
  auto v_hocn = atlas::array::make_view<double, 2>(xb["hocn"]);
  auto v_tocn = atlas::array::make_view<double, 2>(xb["tocn"]);

  // calculate T background error
  atlas::Field tocn_err = innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(levels) |
    atlas::option::name("tocn"));
  bkgErr_.add(tocn_err);
  auto v_tocn_err = atlas::array::make_view<double, 2>(tocn_err);
  {
    v_tocn_err.assign(1.0);
    std::vector<double> dtdz(levels);
    for (size_t i = 0; i < tocn_err.shape(0); i++) {
      if (v_mask(i, 0) == 0.0) continue;
      
      // calculate dt/dz    
      for (size_t z = 1; z < levels-1; z++) {
        dtdz[z] = (v_tocn(i, z+1) - v_tocn(i, z-1)) / 
                  (v_hocn(i,z) + 0.5 * (v_hocn(i, z+1) + v_hocn(i, z-1)));
      }
      dtdz[0] = dtdz[1];
      dtdz[levels-1] = dtdz[levels-2];

      // calculate value as function of dt/dz, efolding scale, and min/max
      for (size_t z = 0; z < levels; z++) {
        v_tocn_err(i, z) = params.tocn.value().min;

        // v_tocn_err(i, z) = abs(params.tocn.value().dz * dtdz[z]);
      }
    }
  }


  // calculate unbalanced S background error

  bkgErr_.add(
    innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(levels) |
    atlas::option::name("socn")));

  // calculate unbalanced SSH background error
  bkgErr_.add(
    innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(1) |
    atlas::option::name("ssh")));

  if (params.output.value() != boost::none) {

    xb.write(*params.output.value());
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
