/*
 * (C) Copyright 2021-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Timer.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/Traits.h"
#include "soca/LinearVariableChange/LinearModel2GeoVaLs/LinearModel2GeoVaLs.h"

namespace soca {

static LinearVariableChangeMaker<LinearModel2GeoVaLs>
       makerLinearVariableChangeModel2GeoVaLs_("LinearModel2GeoVaLs");
static LinearVariableChangeMaker<LinearModel2GeoVaLs>
       makerLinearVariableChangeModel2GeoVaLsDefault_("default");

// -----------------------------------------------------------------------------

LinearModel2GeoVaLs::LinearModel2GeoVaLs(const State & bg, const State &fg,
                                        const Geometry &geom,
                                        const eckit::Configuration &conf)
  : geom_(geom) {
}

// -----------------------------------------------------------------------------

LinearModel2GeoVaLs::~LinearModel2GeoVaLs() {
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiply(const Increment &dxin,
                                         Increment &dxout) const {
  util::Timer timer("soca::LinearModel2GeoVaLs", "multiply");

  const auto & fsetIn = dxin.fieldSet();
  auto & fsetOut = dxout.fieldSet();

  // Identity operator. The only thing this does is check for names being
  // different, and if a 2D field is requested from a 3D field then only the
  // surface is retrieved
  for (auto & fOut : fsetOut) {
    // get the correct input field (mapping variable names where necessary)
    const auto & fOutMeta = geom_.fieldMetadata(fOut.name());
    const auto & fIn = fsetIn[fOutMeta.name];

    // copy the data (turning 3D to 2D if necessary)
    const auto & v_fIn = atlas::array::make_view<double, 2>(fIn);
    auto v_fOut = atlas::array::make_view<double, 2>(fOut);
    ASSERT(fOut.shape(1) <= fIn.shape(1));
    for (atlas::idx_t i = 0; i < fOut.shape(0); i++) {
      for (atlas::idx_t lvl = 0; lvl < fOut.shape(1); lvl++) {
        v_fOut(i, lvl) = v_fIn(i, lvl);
      }
    }
  }
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiplyInverse(const Increment &dxin,
                                                Increment &dxout) const {
  util::Timer timer("soca::LinearModel2GeoVaLs", "multiplyInverse");
  multiply(dxin, dxout);
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiplyAD(const Increment &dxin,
                                           Increment &dxout) const {
  util::Timer timer("soca::LinearModel2GeoVaLs", "multiplyAD");

  const auto & fsetIn = dxin.fieldSet();
  auto & fsetOut = dxout.fieldSet();

  // adjoint. Add fields, with identity transformation, to the output
  for (auto & fIn : fsetIn) {
    // get the correct output field (mapping variable names where necessary)
    const auto & fInMeta = geom_.fieldMetadata(fIn.name());
    auto & fOut = fsetOut[fInMeta.name];

    // add data to the output
    const auto & v_fIn = atlas::array::make_view<double, 2>(fIn);
    auto v_fOut = atlas::array::make_view<double, 2>(fOut);
    ASSERT(fOut.shape(1) >= fIn.shape(1));
    for (atlas::idx_t i = 0; i < fOut.shape(0); i++) {
      for (atlas::idx_t lvl = 0; lvl < fIn.shape(1); lvl++) {
        v_fOut(i, lvl) += v_fIn(i, lvl);
      }
    }
  }
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::multiplyInverseAD(const Increment &dxin,
                                                  Increment &dxout) const {
  util::Timer timer("soca::LinearModel2GeoVaLs", "multiplyInverseAD");
  multiplyAD(dxin, dxout);
}

// -----------------------------------------------------------------------------

void LinearModel2GeoVaLs::print(std::ostream & os) const {
  os << "SOCA linear change variable: LinearModel2GeoVaLs";
}

// -----------------------------------------------------------------------------

}  // namespace soca
