/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "soca/SaberBlocks/ParametricOceanStdDev/ParametricOceanStdDev.h"

namespace soca {

// ----------------------------------------------------------------------------------------

static saber::SaberOuterBlockMaker<ParametricOceanStdDev> makerParametricOceanStdDev_("SOCAParametricOceanStdDev");

// ----------------------------------------------------------------------------------------

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
  util::Timer timer("soca::ParametricOceanStdDev", "ParametricOceanStdDev");
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiply(oops::FieldSet3D &fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::multiply starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "multiply");

  // TODO implement this

  oops::Log::trace() << "ParametricOceanStdDev::multiply done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::multiplyAD starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "multiplyAD");

  // TODO implement this

  oops::Log::trace() << "ParametricOceanStdDev::multiplyAD done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::leftInverseMultiply(oops::FieldSet3D &fset) const {
  oops::Log::trace() << "ParametricOceanStdDev::leftInverseMultiply starting" << std::endl;
  util::Timer timer("soca::ParametricOceanStdDev", "leftInverseMultiply");

  // TODO implement this

  oops::Log::trace() << "ParametricOceanStdDev::leftInverseMultiply done" << std::endl;
}

// ----------------------------------------------------------------------------------------

void ParametricOceanStdDev::print(std::ostream & os) const {
  os << "soca::ParametricOceanStdDev";
}

// ----------------------------------------------------------------------------------------

}  // namespace soca
