/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

// -----------------------------------------------------------------------------------------
namespace oops{

  class FieldSet3D;
}
namespace eckit {
  class Configuration;
}

// -----------------------------------------------------------------------------------------

namespace soca
{

class ParametricOceanStdDevParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ParametricOceanStdDevParameters, saber::SaberBlockParametersBase)
 public:
   oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------------------

class ParametricOceanStdDev : public saber::SaberOuterBlockBase {
 public:
  static const std::string classname() {return "soca::ParametricOceanStdDev";}

  typedef ParametricOceanStdDevParameters Parameters_;

  ParametricOceanStdDev( const oops::GeometryData &,
    const oops::Variables &,
    const eckit::Configuration &,
    const Parameters_ &,
    const oops::FieldSet3D &,
    const oops::FieldSet3D &);
  
  virtual ~ParametricOceanStdDev() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
   void print(std::ostream &) const override;

   const oops::GeometryData & innerGeometryData_;
   oops::Variables innerVars_;
};

    
} // namespace soca
