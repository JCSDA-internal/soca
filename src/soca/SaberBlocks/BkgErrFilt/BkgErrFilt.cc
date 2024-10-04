/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <set>

#include "oops/util/Timer.h"

#include "soca/SaberBlocks/BkgErrFilt/BkgErrFilt.h"

namespace soca
{

// ----------------------------------------------------------------------------------------

static saber::SaberOuterBlockMaker<BkgErrFilt> makerBkgErrFilt_("SOCABkgErrFilt");

// ----------------------------------------------------------------------------------------

BkgErrFilt::BkgErrFilt(
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
  util::Timer timer("soca::BkgErrFilt", "BkgErrFilt");

  // TODO(travis) this block is not working with dual resolution.
  // outerGeometryData is the low resolution grid, but xb and fg are on the
  // high resolution grid. I don't think I have what I need to to interpolate.
  ASSERT(outerGeometryData.fieldSet()[0].shape(0) == xb.fieldSet()[0].shape(0));

  const double MIN_THICKNESS = 1.0e-3;
  const int levels = xb["hocn"].levels();

  mask_ = innerGeometryData_.getField("interp_mask");
  mult3D_ = innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(levels) );
  mult2D_ = innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(1) );
  auto v_mult3D = atlas::array::make_view<double, 2>(mult3D_);
  auto v_mult2D = atlas::array::make_view<double, 2>(mult2D_);
  const auto v_hocn = atlas::array::make_view<double, 2>(xb["hocn"]);
  const auto v_mask = atlas::array::make_view<double, 2>(mask_);

  v_mult3D.assign(0.0);
  v_mult2D.assign(0.0);
  for (size_t i = 0; i < innerGeometryData_.functionSpace().size(); i++) {
    if (v_mask(i, 0) == 0.0) continue;

    // mask out if total depth is less than oceanDepthMin
    double depth = 0.0;
    for (size_t z = 0; z < levels; z++) depth += v_hocn(i, z);
    if (depth <= params.oceanDepthMin) continue;

    // 2D vars (i.e. ssh) get masked and rescaled by a constant
    v_mult2D(i, 0) = params.rescaleBkgerr;

    // 3D vars (i.e. socn, tocn) get masked and rescaled by an exponential decay with depth
    depth = 0.0;
    for (size_t z = 0; z < levels; z++) {
      depth += v_hocn(i, z) / 2.0;  // keep track of layer depth as we go
      if (v_hocn(i, z) > MIN_THICKNESS) {
        v_mult3D(i, z) = params.rescaleBkgerr * exp(-depth/params.efoldZ);
      }
      depth += v_hocn(i, z) / 2.0;  // add the other half of the layer for next loop iteration
    }
  }
}

// ----------------------------------------------------------------------------------------

void BkgErrFilt::multiply(oops::FieldSet3D & fset) const {
  util::Timer timer("soca::BkgErrFilt", "multiply");

  const std::set<std::string> FIELDS3D{"sea_water_potential_temperature", "socn"};
  const std::set<std::string> FIELDS2D{"sea_surface_height_above_geoid"};

  const auto v_mult3D = atlas::array::make_view<double, 2>(mult3D_);
  const auto v_mult2D = atlas::array::make_view<double, 2>(mult2D_);
  const auto v_mask   = atlas::array::make_view<double, 2>(mask_);

  for (atlas::Field field : fset) {
    const auto levels = field.shape(1);
    auto v_field = atlas::array::make_view<double, 2>(field);
    const bool inFields3D = FIELDS3D.find(field.name()) != FIELDS3D.end();
    const bool inFields2D = FIELDS2D.find(field.name()) != FIELDS2D.end();

    for (size_t i = 0; i < field.shape(0); i++) {
      // Note, the IF statements are outside the loop over z to avoid the
      // overhead of checking the same condition for every z level.
      if (inFields2D) {
        // 2D vars (i.e. ssh) get masked and rescaled by a constant
        v_field(i, 0) *= v_mult2D(i, 0);
      } else if (inFields3D) {
        // 3D vars (i.e. socn, tocn) get masked and rescaled by an exponential decay with depth
        for (size_t z = 0; z < levels; z++) {
          v_field(i, z) *= v_mult3D(i, z);
        }
      } else {
        // all other vars get masked
        for (size_t z = 0; z < levels; z++) {
          v_field(i, z) *= v_mask(i, 0);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------------------

void BkgErrFilt::multiplyAD(oops::FieldSet3D & fset) const {
  multiply(fset);
}

// ----------------------------------------------------------------------------------------

void BkgErrFilt::print(std::ostream & os) const {
  os << classname();
}

// ----------------------------------------------------------------------------------------

}  // namespace soca
