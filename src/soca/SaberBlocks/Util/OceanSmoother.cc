/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/SaberBlocks/Util/OceanSmoother.h"
#include "soca/Utils/Diffusion.h"

#include "atlas/field.h"

#include "oops/util/Logger.h"

namespace soca
{

// --------------------------------------------------------------------------------------

OceanSmoother::OceanSmoother(const oops::GeometryData & geom, const Parameters_ & params, const int levels)
 :  geom_(geom), diffusion_(geom)
{
  atlas::FieldSet scales;

  // helper function to mask the scales
  // TODO include a minimum layer thickness so that bottom layers are ignored
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
          v_hzScales(i, 0) += v_rossby(i,0) * rossbyMult;
      }
    }

    // impose min based on grid size
    const double minGridMult = hzParam.minGridMult.value();
    oops::Log::info() << "  minimum grid-based multiplier: " << minGridMult << std::endl;
    if (minGridMult > 0.0) {
      const auto &v_area = atlas::array::make_view<double, 2>(geom_.getField("area"));
      for (size_t i = 0; i < hzScales.shape(0); i++) {
          v_hzScales(i, 0) = std::max(v_hzScales(i, 0), std::sqrt(v_area(i,0))*minGridMult);
      }
    }

    // impose global min/max
    const double minVal = hzParam.min.value();
    const double maxVal = hzParam.max.value();
    oops::Log::info() << "  global minimum: " << minVal << std::endl;
    oops::Log::info() << "  global maximum: " << maxVal << std::endl;
    for (size_t i = 0; i < hzScales.shape(0); i++) {
        v_hzScales(i, 0) = std::clamp(v_hzScales(i, 0), minVal, maxVal);
    }

    // optional mask
    oops::Log::info() << "  using mask: "
                      << (params.mask.value() == "" ? "NONE" : params.mask.value() ) << std::endl;
    maskField(hzScales);
  }

  // calculate vertical scales
  if (params.vertical.value() != boost::none) {
    const OceanSmootherParameters::Vertical & vtParam = *params.vertical.value();
    atlas::Field vtScales = geom.functionSpace().createField<double>(
      atlas::option::name("vtScales") | atlas::option::levels(levels));
    scales.add(vtScales);
    auto v_vtScales = atlas::array::make_view<double, 2>(vtScales);

    oops::Log::info() << "Ocean Smoother : calculating horizontal scales: " << std::endl;

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

    // optional mask
    oops::Log::info() << "  using mask: "
                      << (params.mask.value() == "" ? "NONE" : params.mask.value() ) << std::endl;
    maskField(vtScales);

    // TODO zero out thin layers at the bottom

    // TODO smooth by the horizontal scales
  }

  // all done, set the scales of the diffusion smoother
  diffusion_.setScales(scales);
}

// --------------------------------------------------------------------------------------

void OceanSmoother::multiply(atlas::FieldSet & fset) const {
  diffusion_.multiply(fset, Diffusion::HZVT_2D_1D, true);
}

// --------------------------------------------------------------------------------------

void OceanSmoother::multiply(atlas::Field & field) const {
  diffusion_.multiply(field, Diffusion::HZVT_2D_1D, true);
}

}