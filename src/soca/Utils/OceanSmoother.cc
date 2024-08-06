/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "atlas/field.h"

#include "oops/base/GeometryData.h"
#include "oops/generic/Diffusion.h"
#include "oops/util/FieldSetHelpers.h"

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
  atlas::FieldSet scales;  // The final fields sent to the diffusion operator
  atlas::FieldSet diags;  // Diagnostic fields that are optionally output
  const auto & v_ghost = atlas::array::make_view<int, 1>(geom_.functionSpace().ghost());

  // Create the 3D mask field
  {
    mask_.reset(new atlas::Field(geom.functionSpace().createField<double>(
      atlas::option::name("mask") | atlas::option::levels(levels))));
    auto v_mask = atlas::array::make_view<double, 2>(*mask_);
    v_mask.assign(1.0);
    diags.add(*mask_);

    // apply 2D horizontal mask
    const auto & maskParams = params.mask.value();
    if (maskParams.maskVariable.value() != "") {
      const auto & v_hzMask = atlas::array::make_view<double, 2>(
        geom_.getField(maskParams.maskVariable.value()));
      for (size_t i = 0; i < mask_->shape(0); i++) {
        if (v_ghost(i)) continue;
        for (size_t level = 0; level < mask_->shape(1); level++) {
          v_mask(i, level) = v_hzMask(i, 0) ? 1.0 : 0.0;
        }
      }
      mask_->set_dirty();
    }

    // mask out vertical layers that are too thin
    const auto minThickness = maskParams.minThickness.value();
    if (minThickness > 0.0) {
      const auto & v_thickness = atlas::array::make_view<double, 2>(
        bkg.field(params.thicknessVariable.value()));
      for (size_t i = 0; i < mask_->shape(0); i++) {
        if (v_ghost(i)) continue;
        for (size_t level = 0; level < mask_->shape(1); level++) {
          if (v_thickness(i, level) < minThickness) {
            v_mask(i, level) = 0.0;
          }
        }
      }
      mask_->set_dirty();
    }
  }

  // helper function to mask the scales later
  auto maskField = [this, &bkg, &params, &v_ghost](atlas::Field &field) {
    if (mask_) {
      const auto & v_mask = atlas::array::make_view<double, 2>(*mask_);
      auto v_field = atlas::array::make_view<double, 2>(field);
      for (size_t i = 0; i < field.shape(0); i++) {
        if (v_ghost(i)) continue;
        for (size_t level = 0; level < field.shape(1); level++) {
          v_field(i, level) *= v_mask(i, level);
        }
      }
      field.set_dirty();
    }
  };

  // calculate horizontal scales (in meters)
  if (params.horizontal.value() != boost::none) {
    const OceanSmoother::Parameters::Horizontal & hzParam = *params.horizontal.value();
    atlas::Field hzScales = geom.functionSpace().createField<double>(
      atlas::option::name("hzScales") | atlas::option::levels(levels));
    hzScales->set_dirty();
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
        if (v_ghost(i)) continue;
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
        if (v_ghost(i)) continue;
        for (size_t level = 0; level < levels; level++) {
          v_hzScales(i, level) = std::max(
            v_hzScales(i, level),
            std::sqrt(v_area(i, 0))*minGridMult);
        }
      }
    }

    // impose global min/max
    const double minVal = hzParam.min.value();
    const double maxVal = hzParam.max.value();
    oops::Log::info() << "  global minimum: " << minVal << std::endl;
    oops::Log::info() << "  global maximum: " << maxVal << std::endl;
    for (size_t i = 0; i < hzScales.shape(0); i++) {
      if (v_ghost(i)) continue;
      for (size_t level = 0; level < levels; level++) {
        v_hzScales(i, level) = std::clamp(v_hzScales(i, level), minVal, maxVal);
      }
    }

    // optional mask
    maskField(hzScales);

    // done with horizontal scales
    scales.add(hzScales);
    diags.add(hzScales);
  }

  // calculate vertical scales (number of levels)
  if (params.vertical.value() != boost::none) {
    const OceanSmoother::Parameters::Vertical & vtParam = *params.vertical.value();
    atlas::Field vtScales = geom.functionSpace().createField<double>(
      atlas::option::name("vtScales") | atlas::option::levels(levels));
    auto v_vtScales = atlas::array::make_view<double, 2>(vtScales);
    vtScales.set_dirty();
    v_vtScales.assign(0.0);

    oops::Log::info() << "Ocean Smoother : calculating vertical scales: " << std::endl;

    // use the MLD at the surface do an optional MLD based smoothing
    if (vtParam.useMld.value()) {
      oops::Log::info() << "  using MLD for vertical scales" << std::endl;
      // get MLD
      const std::string & mldVariable = vtParam.mldVariable.value();
      auto mld = bkg.field(mldVariable).clone();  // make a copy, because we'll modify it
      diags.add(mld);
      auto v_mld = atlas::array::make_view<double, 2>(mld);

      // smooth the MLD field by the horizontal scales, if requested
      if (vtParam.mldSmooth.value()) {
        oops::Log::info() << "  smoothing MLD field" << std::endl;
        atlas::FieldSet fset;
        fset.add(mld);
        oops::Diffusion hzSmoother(*diffusion_);
        hzSmoother.setParameters(scales);
        hzSmoother.multiply(fset, oops::Diffusion::Mode::HorizontalOnly);
      }

      // cap max MLD
      const double mldMax = vtParam.mldMax.value();
      oops::Log::info() << "  maximum MLD value: " << mldMax << std::endl;
      for (size_t i = 0; i < mld.shape(0); i++) {
        if (v_ghost(i)) continue;
        v_mld(i, 0) = std::min(v_mld(i, 0), mldMax);
      }

      // calculate the fractional number of levels in the MLD
      const auto & v_mask = atlas::array::make_view<double, 2>(*mask_);
      const auto & v_hocn = atlas::array::make_view<double, 2>(
        bkg.field(params.thicknessVariable.value()));
      auto mldLevels = mld.clone();
      mldLevels.rename("mldLevels");
      diags.add(mldLevels);
      auto v_mldLevels = atlas::array::make_view<double, 2>(mldLevels);
      v_mldLevels.assign(0.0);
      for (size_t i = 0; i < mldLevels.shape(0); i++) {
        if (v_ghost(i)) continue;
        if (v_mask(i, 0) == 0.0) continue;  // skip land points

        // find the last level in the MLD, the depth at that level, and the depth of the next level
        double depth = 0.0;
        double depth_p1 = 0.0;
        double mldLevel = levels - 2;
        for (size_t level = 0; level < levels; level++) {
          depth_p1 += v_hocn(i, level)/2.0;
          if (v_mask(i, level) == 0.0 || depth_p1 > v_mld(i, 0)) {
            mldLevel = 1.0*level - 1;
            break;
          }
          depth = depth_p1;
          depth_p1 += v_hocn(i, level)/2.0;
        }

        // what percentage of next level to add. Note, we clamp to [0,1] because
        // the MLD could be deeper than the last level due to the earlier
        // smoothing of MLD
        double frac = 1.0;
        if (depth_p1 > depth) {
          frac = std::clamp((v_mld(i, 0) - depth) /(depth_p1 - depth), 0.0, 1.0);
        }

        // final
        v_mldLevels(i, 0) = mldLevel + frac;
      }

      // calculate the vertical scales based on the MLD
      for (size_t i = 0; i < vtScales.shape(0); i++) {
        if (v_ghost(i)) continue;
        for (size_t level = 0; level < levels; level++) {
          v_vtScales(i, level) = v_mldLevels(i, 0)-level;
        }
      }
    }

    // impose global min/max
    const double minVal = vtParam.min.value();
    const double maxVal = vtParam.max.value();
    oops::Log::info() << "  global minimum: " << minVal << std::endl;
    oops::Log::info() << "  global maximum: " << maxVal << std::endl;
    ASSERT(minVal <= maxVal);
    for (size_t i = 0; i < vtScales.shape(0); i++) {
      if (v_ghost(i)) continue;
      for (size_t level = 0; level < levels; level++) {
        v_vtScales(i, level) = std::clamp(v_vtScales(i, level), minVal, maxVal);
      }
    }

    // optional mask
    maskField(vtScales);

    // done with vertical scales
    scales.add(vtScales);
    diags.add(vtScales);
  }

  // Done, set the parameters and initialize the diffusion class
  diffusion_->setParameters(scales);

  // save diagnostics, if requested
  if (params.saveDiags.value() != boost::none) {
    oops::Log::info() << "Ocean Smoother : saving diagnostics to file " << std::endl;
    util::writeFieldSet(geom_.comm(), *params.saveDiags.value(), diags);
  }

  // keep the 3D mask around ONLY if we need to set masked values to 0 in
  // multiply, otherwise delete it to save memory
  if (!params.mask.value().setMaskedToZero.value()) {
    mask_.reset();
  }
}

// -----------------------------------------------------------------------------

void OceanSmoother::multiply(atlas::FieldSet & fields) {
  // diffusion based 3D smoothing
  diffusion_->multiply(fields);

  // mask out the fields if requested
  if (mask_) {
    const auto & v_ghost = atlas::array::make_view<int, 1>(geom_.functionSpace().ghost());
    for (size_t f = 0; f < fields.size(); f++) {
      auto & field = fields[f];
      auto v_field = atlas::array::make_view<double, 2>(field);
      auto v_mask = atlas::array::make_view<double, 2>(*mask_);
      for (size_t i = 0; i < field.shape(0); i++) {
        if (v_ghost(i)) continue;
        for (size_t level = 0; level < field.shape(1); level++) {
          v_field(i, level) *= v_mask(i, level);
        }
      }
      field.set_dirty();
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

}  // namespace soca
