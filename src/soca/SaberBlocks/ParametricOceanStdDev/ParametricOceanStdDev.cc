/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/SaberBlocks/ParametricOceanStdDev/ParametricOceanStdDev.h"

#include "saber/blocks/SaberOuterBlockBase.h"

namespace soca {

static saber::SaberOuterBlockMaker<ParametricOceanStdDev> makerParametricOceanStdDev_("ParametricOceanStdDev");

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
  ASSERT(1==0);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiply(oops::FieldSet3D &) const
{
  ASSERT(1==2);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiplyAD(oops::FieldSet3D &) const
{
  ASSERT(1==2);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::leftInverseMultiply(oops::FieldSet3D &) const
{
  ASSERT(1==2);
}

// ------------------------------------------------------------------------------------------------

void ParametricOceanStdDev::print(std::ostream & os) const {
  os << classname();
}

// ------------------------------------------------------------------------------------------------
}