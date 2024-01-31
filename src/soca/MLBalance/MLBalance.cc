#include <typeinfo>

#include "atlas/array.h"

#include "soca/MLBalance/MLBalance.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"

namespace soca {

// --------------------------------------------------------------------------------------

static saber::SaberOuterBlockMaker<MLBalance> makerMLBalance_("MLBalance");

// --------------------------------------------------------------------------------------

MLBalance::MLBalance(
    const oops::GeometryData & geometryData,
    const oops::Variables & vars,
    const eckit::Configuration & mlbConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberOuterBlockBase(params, xb.validTime()),
    params_(params), innerVars_(vars), innerGeometryData_(geometryData), kct_(0.1)
{
  // Get active variables
  vars_ = params.activeVars.value().get_value_or(vars);

  // Do something with the ML Balances
  oops::Log::info() << params_.mlbalances.value() << std::endl;
}

// --------------------------------------------------------------------------------------

MLBalance::~MLBalance() {}

// --------------------------------------------------------------------------------------

void MLBalance::multiply(oops::FieldSet3D & fset) const {
  // get the indices corresponding to cicen and tocn
  // TODO (G): There's probably a better way to do this?
  int index = 0;
  int dcIndex;
  int dtIndex;
  for (auto & field : fset) {
    if (field.name() == "cicen") {
      dcIndex = index;
    }
    if (field.name() == "tocn") {
      dtIndex = index;
    }
    index += 1;
  }

  // Apply forward balance
  // dc +=  kct * dsst
  auto dcView = atlas::array::make_view<double, 2>(fset[dcIndex]);
  auto dtView = atlas::array::make_view<double, 2>(fset[dtIndex]);
  for (int jnode = 0; jnode < fset[dcIndex].shape(0); ++jnode) {
      dcView(jnode, 0) += kct_ * dtView(jnode, 0);
  }
}

// --------------------------------------------------------------------------------------

void MLBalance::multiplyAD(oops::FieldSet3D & fset) const {
  // get the indices corresponding to cicen and tocn
  // TODO (G): There's probably a better way to do this?
  int index = 0;
  int dcIndex;
  int dtIndex;
  for (auto & field : fset) {
    if (field.name() == "cicen") {
      dcIndex = index;
    }
    if (field.name() == "tocn") {
      dtIndex = index;
    }
    index += 1;
  }

  // Apply adjoint balance
  // tocn +=  kct * cicen
  auto dcView = atlas::array::make_view<double, 2>(fset[dcIndex]);
  auto dtView = atlas::array::make_view<double, 2>(fset[dtIndex]);
  for (int jnode = 0; jnode < fset[dcIndex].shape(0); ++jnode) {
    dtView(jnode, 0) += kct_ * dcView(jnode, 0);
  }
}

// --------------------------------------------------------------------------------------

void MLBalance::print(std::ostream &) const {
}

// --------------------------------------------------------------------------------------

}  // namespace soca
