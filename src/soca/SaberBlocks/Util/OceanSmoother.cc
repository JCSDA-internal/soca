/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/SaberBlocks/Util/OceanSmoother.h"
#include "soca/Utils/Diffusion.h"

#include "atlas/field.h"

namespace soca
{

// --------------------------------------------------------------------------------------

OceanSmoother::OceanSmoother(const oops::GeometryData & geom, const Parameters_ & params)
 :  geom_(geom), diffusion_(geom)
{
  atlas::FieldSet scales;

  // helper function to mask the scales
  auto maskField=[this, &params](atlas::Field &field) {
    if (params.mask.value() == "") return;

    const auto & v_mask = atlas::array::make_view<double, 2>(geom_.getField(params.mask.value()));
    auto v_field = atlas::array::make_view<double, 2>(field);
    for (size_t i = 0; i < field.shape(0); i++) {
      if (!v_mask(i,0)) {
        for(size_t level = 0; level < field.shape(1); level++) {
          v_field(i, level) = 0.0;
        }
      }
    }
  };

  // calculate horizontal scales
  if (params.horizontal.value() != boost::none) {
    const OceanSmootherParameters::Horizontal & hzParam = *params.horizontal.value();
    atlas::Field hzScales = geom.functionSpace().createField<double>(
      atlas::option::name("hzScales") | atlas::option::levels(1));
    scales.add(hzScales);
    auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);

    // base value
    v_hzScales.assign(hzParam.base.value());

    // + rossby radius based value
    if (hzParam.rossbyMult.value() > 0.0) {
      const std::string & rossbyVariable = hzParam.rossbyVariable.value();
      const auto &v_rossby = atlas::array::make_view<double, 2>(geom_.getField(rossbyVariable));
      for (size_t i = 0; i < hzScales.shape(0); i++) {
          v_hzScales(i, 0) += v_rossby(i,0) * hzParam.rossbyMult.value();
      }
    }

    // impose min based on grid size
    if (hzParam.minGridMult.value() > 0.0) {
      const auto &v_area = atlas::array::make_view<double, 2>(geom_.getField("area"));
      for (size_t i = 0; i < hzScales.shape(0); i++) {
          v_hzScales(i, 0) = std::max(v_hzScales(i, 0), std::sqrt(v_area(i,0))*hzParam.minGridMult.value());
      }
    }

    // impose global min/max
    const double minVal = hzParam.min.value();
    const double maxVal = hzParam.max.value();
    for (size_t i = 0; i < hzScales.shape(0); i++) {
        v_hzScales(i, 0) = std::clamp(v_hzScales(i, 0), minVal, maxVal);
    }

    // optional mask
    maskField(hzScales);
  }

  // calculate vertical scales
  if (params.vertical.value() != boost::none) {
    const OceanSmootherParameters::Vertical & vtParam = *params.vertical.value();
    atlas::Field vtScales = geom.functionSpace().createField<double>(
      atlas::option::name("vtScales") | atlas::option::levels(vtParam.levels.value()));
    scales.add(vtScales);
    auto v_vtScales = atlas::array::make_view<double, 2>(vtScales);

    v_vtScales.assign(vtParam.base.value());


    // TODO set the rest

    maskField(vtScales);
  }

  // all done, set the scales of the diffusion smoother
  diffusion_.setScales(scales);
}

// --------------------------------------------------------------------------------------

void OceanSmoother::multiply(atlas::FieldSet & fset) const {
  diffusion_.multiply(fset);
}

// --------------------------------------------------------------------------------------

void OceanSmoother::multiply(atlas::Field & field) const {
  diffusion_.multiply(field);
}

}