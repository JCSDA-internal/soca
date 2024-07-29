/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/Diffusion.h"

#include "soca/Utils/OceanSmoother.h"

namespace soca {

// -----------------------------------------------------------------------------

OceanSmoother::OceanSmoother(
    const oops::GeometryData & geom,
    const Parameters & params,
    const int levels,
    const atlas::FieldSet & bkg)
  : geom_(geom),
    diffusion_(new oops::Diffusion(geom_))
{
  atlas::FieldSet scales;

  // Create the 3D mask field
  {
    mask_.reset(new atlas::Field(geom.functionSpace().createField<double>(
      atlas::option::name("mask") | atlas::option::levels(levels))));
    auto v_mask = atlas::array::make_view<double, 2>(*mask_);
    v_mask.assign(1.0); // true where ocean is valid

    const auto maskParams = params.mask.value();

    // 2D horizontal mask
    if (maskParams.maskVariable.value() != "") {
      const auto & v_hzMask = atlas::array::make_view<double, 2>(geom_.getField(maskParams.maskVariable.value()));
      for (size_t i = 0; i < mask_->shape(0); i++) {
        for(size_t level = 0; level < mask_->shape(1); level++) {
          v_mask(i, level) = v_hzMask(i,0) ? 1.0 : 0.0;
        }
      }
    }

    // mask out layers that are too thin
    const auto minThickness = maskParams.minThickness.value();
    if (minThickness > 0.0) {
      const auto & v_thickness = atlas::array::make_view<double, 2>(bkg.field(maskParams.thicknessVariable.value()));
      for (size_t i = 0; i < mask_->shape(0); i++) {
        for(size_t level = 0; level < mask_->shape(1); level++) {
          if (v_thickness(i, level) < minThickness) {
            v_mask(i, level) = 0.0;
          }
        }
      }
    }
  }

  // helper function to mask the scales later
  auto maskField=[this, &bkg, &params](atlas::Field &field) {
    if (mask_) {
      const auto & v_mask = atlas::array::make_view<double, 2>(*mask_);
      auto v_field = atlas::array::make_view<double, 2>(field);
      for (size_t i = 0; i < field.shape(0); i++) {
        for(size_t level = 0; level < field.shape(1); level++) {
          v_field(i, level) *= v_mask(i, level);
        }
      }
    }
  };

  // calculate horizontal scales
  if (params.horizontal.value() != boost::none) {
    const OceanSmoother::Parameters::Horizontal & hzParam = *params.horizontal.value();
    atlas::Field hzScales = geom.functionSpace().createField<double>(
      atlas::option::name("hzScales") | atlas::option::levels(levels));
    scales.add(hzScales);
    auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);

    oops::Log::info() << "Ocean Smoother : calculating horizontal scales: " << std::endl;

    // base value
    const double base = hzParam.base.value();
    oops::Log::info() << "  base value: " << base << std::endl;
    v_hzScales.assign(base);

    // + rossby radius based value
    const double rossbyMult = hzParam.rossbyMult.value();
    oops::Log::info() << "  rossby radius multiplier: " << rossbyMult << std::endl;
    if (rossbyMult > 0.0) {
      const std::string & rossbyVariable = hzParam.rossbyVariable.value();
      const auto &v_rossby = atlas::array::make_view<double, 2>(geom_.getField(rossbyVariable));
      for (size_t i = 0; i < hzScales.shape(0); i++) {
        for (size_t level = 0; level < levels; level++) {
          v_hzScales(i, level) += v_rossby(i, level) * rossbyMult;
        }
      }
    }

    // impose min based on grid size
    const double minGridMult = hzParam.minGridMult.value();
    oops::Log::info() << "  minimum grid-based multiplier: " << minGridMult << std::endl;
    if (minGridMult > 0.0) {
      const auto &v_area = atlas::array::make_view<double, 2>(geom_.getField("area"));
      for (size_t i = 0; i < hzScales.shape(0); i++) {
        for (size_t level = 0; level < levels; level++) {
          v_hzScales(i, level) = std::max(v_hzScales(i, level), std::sqrt(v_area(i,0))*minGridMult);
        }
      }
    }

    // impose global min/max
    const double minVal = hzParam.min.value();
    const double maxVal = hzParam.max.value();
    oops::Log::info() << "  global minimum: " << minVal << std::endl;
    oops::Log::info() << "  global maximum: " << maxVal << std::endl;
    for (size_t i = 0; i < hzScales.shape(0); i++) {
      for (size_t level = 0; level < levels; level++) {
        v_hzScales(i, level) = std::clamp(v_hzScales(i, level), minVal, maxVal);
      }
    }

    // optional mask
    maskField(hzScales);
  }

  // calculate vertical scales
  if (params.vertical.value() != boost::none) {
    const OceanSmoother::Parameters::Vertical & vtParam = *params.vertical.value();
    atlas::Field vtScales = geom.functionSpace().createField<double>(
      atlas::option::name("vtScales") | atlas::option::levels(levels));
    scales.add(vtScales);
    auto v_vtScales = atlas::array::make_view<double, 2>(vtScales);

    oops::Log::info() << "Ocean Smoother : calculating vertical scales: " << std::endl;

    // base value
    const double base = vtParam.base.value();
    oops::Log::info() << "  base value: " << base << std::endl;
    v_vtScales.assign(base);

    // TODO, do an optional MLD based smoothing
    // // impose global min/max
    // // impose global min/max
    // const double minVal = vtParam.min.value();
    // const double maxVal = vtParam.max.value();
    // oops::Log::info() << "  global minimum: " << minVal << std::endl;
    // oops::Log::info() << "  global maximum: " << maxVal << std::endl;
    // for (size_t i = 0; i < vtScales.shape(0); i++) {
    //   for (size_t lvl = 0; lvl < vtScales.shape(1); lvl++) {
    //     v_vtScales(i, lvl) = std::clamp(v_vtScales(i, lvl), minVal, maxVal);
    //   }
    // }

    // TODO optional smoothing of vt scales by the horizontal scales

    // optional mask
    maskField(vtScales);
  }

  // Done, set the parameters and initialize the diffusion class
  diffusion_->setParameters(scales);

  // keep the 3D mask around ONLY if we need to set masked values to 0 in
  // multiply
  if (!params.mask.value().setMaskedToZero.value()) {
    mask_.reset();
  }
}

// -----------------------------------------------------------------------------

void OceanSmoother::multiply(atlas::FieldSet & fields) {
  // diffusion based 3D smoothing
  diffusion_->multiply(fields);

  // mask out the fields if needed
  if (mask_) {
    for (size_t f = 0; f < fields.size(); f++) {
      auto & field = fields[f];
      auto v_field = atlas::array::make_view<double, 2>(field);
      auto v_mask = atlas::array::make_view<double, 2>(*mask_);
      for (size_t i = 0; i < field.shape(0); i++) {
        for(size_t level = 0; level < field.shape(1); level++) {
          v_field(i, level) *= v_mask(i, level);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void OceanSmoother::multiply(atlas::Field & field) {
  atlas::FieldSet fset;
  fset.add(field);
  multiply(fset);
}

// -----------------------------------------------------------------------------

}