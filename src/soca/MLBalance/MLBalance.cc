/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <typeinfo>

#include "atlas/array.h"

#include "oops/util/FieldSetHelpers.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/MLBalance/MLBalance.h"
#include "soca/MLBalance/MLJac.h"

namespace soca {

// --------------------------------------------------------------------------------------

static saber::SaberOuterBlockMaker<MLBalance> makerMLBalance_("MLBalance");

// --------------------------------------------------------------------------------------

MLBalance::MLBalance(
    const oops::GeometryData & outerGeometryData,
    const oops::Variables & outerVars,
    const eckit::Configuration & mlbConf,
    const Parameters_ & params,
    const oops::FieldSet3D & xb,
    const oops::FieldSet3D & fg)
  : saber::SaberOuterBlockBase(params, xb.validTime()),
    params_(params),
    activeVars_(params.mandatoryActiveVars()),
    innerVars_(outerVars),
    innerGeometryData_(outerGeometryData),
    jac_(nullptr)
{
  // Initialize the Jacobian variables
  std::vector<std::string>
    jacStr{"ds/dt",
           "dssh/dt", "dssh/ds",
           "dc/dsst", "dc/dsss", "dc/dhi", "dc/dhs"};
  int nz(xb["sea_water_potential_temperature"].shape(1));
  const int SINGLE_LEVEL = 1;
  std::vector<int> jacLevels{nz,
                             nz, nz,
                             SINGLE_LEVEL, SINGLE_LEVEL, SINGLE_LEVEL, SINGLE_LEVEL};
  oops::Variables jacVars(jacStr);
  for (size_t i = 0; i < jacStr.size(); ++i) {
    jacVars[jacStr[i]].setLevels(jacLevels[i]);
  }

  // Initialize the Jacobian
  const double INITIAL_VALUE = 0.0;
  jac_ = util::createFieldSet(outerGeometryData.functionSpace(), jacVars, INITIAL_VALUE);

  // Initialize Jacobian
  setupJac(xb, outerGeometryData.comm(), mlbConf);
}

// --------------------------------------------------------------------------------------

MLBalance::~MLBalance() {}

// --------------------------------------------------------------------------------------

void MLBalance::setupJac(const oops::FieldSet3D & xb,
                         const eckit::mpi::Comm & comm,
                         const eckit::Configuration & config) {
  // Create a map of configurations
  eckit::LocalConfiguration mlConf = params_.mlbalances.value();
  MLJac mlJac(mlConf, xb, jac_, innerGeometryData_, comm);
}

// --------------------------------------------------------------------------------------

void MLBalance::multiply(oops::FieldSet3D & fset) const {
  //    Pot. Temp      Salt      ssh     h ice  h snow  ice concentration
  // K= [ I             0         0        0      0      0 ]   Pot. Temp
  //    [ ds/dt         I         0        0      0      0 ]   Salt
  //    [ dssh/dt       dssh/ds   I        0      0      0 ]   ssh
  //    [ 0             0         0        I      0      0 ]   h ice
  //    [ 0             0         0        0      I      0 ]   h snow
  //    [ dc/dsst       dc/dsss   0      dc/dhi dc/dhs   I ]   ice concentration

  // Increment fields
  auto dc = atlas::array::make_view<double, 2>(fset["sea_ice_area_fraction"]);
  auto dhi = atlas::array::make_view<double, 2>(fset["sea_ice_thickness"]);
  auto dhs = atlas::array::make_view<double, 2>(fset["sea_ice_snow_thickness"]);
  auto dt = atlas::array::make_view<double, 2>(fset["sea_water_potential_temperature"]);
  auto ds = atlas::array::make_view<double, 2>(fset["sea_water_salinity"]);
  auto dssh = atlas::array::make_view<double, 2>(fset["sea_surface_height_above_geoid"]);

  // Jacobian fields
  auto dsdt = atlas::array::make_view<double, 2>(jac_["ds/dt"]);
  auto dsshdt = atlas::array::make_view<double, 2>(jac_["dssh/dt"]);
  auto dsshds = atlas::array::make_view<double, 2>(jac_["dssh/ds"]);
  auto dcdsst = atlas::array::make_view<double, 2>(jac_["dc/dsst"]);
  auto dcdsss = atlas::array::make_view<double, 2>(jac_["dc/dsss"]);
  auto dcdhi = atlas::array::make_view<double, 2>(jac_["dc/dhi"]);
  auto dcdhs = atlas::array::make_view<double, 2>(jac_["dc/dhs"]);

  const int nnodes = fset["sea_water_potential_temperature"].shape(0);
  const int nlevs = fset["sea_water_potential_temperature"].shape(1);
  for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
    // Deep copy of some of the input increments
    auto dsshi = dssh(jnode, 0);
    std::vector<double> dsi(ds.shape(1));
    for (size_t j = 0; j < ds.shape(1); ++j) {
      dsi[j] = ds(jnode, j);
    }
    for ( atlas::idx_t jlevel = 0; jlevel < nlevs; ++jlevel ) {
        ds(jnode, jlevel) += dsdt(jnode, jlevel) * dt(jnode, jlevel);
        dssh(jnode, 0) += dsshdt(jnode, jlevel) * dt(jnode, jlevel) +
                          dsshds(jnode, jlevel) * dsi[jlevel];
    }
    dc(jnode, 0) += dcdsst(jnode, 0) * dt(jnode, 0) +
                    dcdsss(jnode, 0) * dsi[0] +
                    dcdhi(jnode, 0) * dhi(jnode, 0) +
                    dcdhs(jnode, 0)* dhs(jnode, 0);
  }
}

// --------------------------------------------------------------------------------------

void MLBalance::multiplyAD(oops::FieldSet3D & fset) const {
  //    Pot. Temp      Salt   ssh     h ice  h snow  ice concentration
  // K^T  = [ I         ds/dt  dssh/dt  0      0      dc/dsst ]   Pot. Temp
  //        [ 0         I      dssh/ds  0      0      dc/dsss ]   Salt
  //        [ 0         0      I        0      0      0       ]   ssh
  //        [ 0         0      0        I      0      dc/dhi  ]   h ice
  //        [ 0         0      0        0      I      dc/dhs  ]   h snow
  //        [ 0         0      0        0      0      I       ]   ice concentration

  // Increment fields
  auto dc = atlas::array::make_view<double, 2>(fset["sea_ice_area_fraction"]);
  auto dhi = atlas::array::make_view<double, 2>(fset["sea_ice_thickness"]);
  auto dhs = atlas::array::make_view<double, 2>(fset["sea_ice_snow_thickness"]);
  auto dt = atlas::array::make_view<double, 2>(fset["sea_water_potential_temperature"]);
  auto ds = atlas::array::make_view<double, 2>(fset["sea_water_salinity"]);
  auto dssh = atlas::array::make_view<double, 2>(fset["sea_surface_height_above_geoid"]);

  // Jacobian fields
  auto dsdt = atlas::array::make_view<double, 2>(jac_["ds/dt"]);
  auto dsshdt = atlas::array::make_view<double, 2>(jac_["dssh/dt"]);
  auto dsshds = atlas::array::make_view<double, 2>(jac_["dssh/ds"]);
  auto dcdsst = atlas::array::make_view<double, 2>(jac_["dc/dsst"]);
  auto dcdsss = atlas::array::make_view<double, 2>(jac_["dc/dsss"]);
  auto dcdhi = atlas::array::make_view<double, 2>(jac_["dc/dhi"]);
  auto dcdhs = atlas::array::make_view<double, 2>(jac_["dc/dhs"]);

  const int nnodes = fset["sea_water_potential_temperature"].shape(0);
  for (atlas::idx_t jnode = 0; jnode < nnodes; ++jnode) {
    auto dci = dc(jnode, 0);
    auto dsshi = dssh(jnode, 0);
    // Deep copy of some of the input increments
    std::vector<double> dsi(ds.shape(1));
    for (size_t j = 0; j < ds.shape(1); ++j) {
      dsi[j] = ds(jnode, j);
    }

    dt(jnode, 0) += dcdsst(jnode, 0) * dci;
    ds(jnode, 0) += dcdsss(jnode, 0) * dci;
    dhi(jnode, 0) += dcdhi(jnode, 0) * dci;
    dhs(jnode, 0) += dcdhs(jnode, 0) * dci;
    const int nlevs = fset["sea_water_potential_temperature"].shape(1);
    for (atlas::idx_t jlevel = 0; jlevel < nlevs; ++jlevel) {
      dt(jnode, jlevel) += dsdt(jnode, jlevel) * dsi[jlevel] +
                           dsshdt(jnode, jlevel) * dsshi;
      ds(jnode, jlevel) += dsshds(jnode, jlevel) * dsshi;
    }
  }
}

// --------------------------------------------------------------------------------------

void MLBalance::print(std::ostream &) const {}

// --------------------------------------------------------------------------------------

}  // namespace soca
