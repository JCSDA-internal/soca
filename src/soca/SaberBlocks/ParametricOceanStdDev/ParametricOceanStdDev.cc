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

  // calculate T background error
  bkgErr_.add(
    innerGeometryData_.functionSpace().createField<double>(
    atlas::option::levels(levels) |
    atlas::option::name("tocn")));

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

}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiply(oops::FieldSet3D &) const
{
  ASSERT(1 == 1);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiplyAD(oops::FieldSet3D &) const
{
  ASSERT(1 == 2);
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
