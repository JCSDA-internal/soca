/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

// ----------------------------------------------------------------------------------------


namespace soca {

// ----------------------------------------------------------------------------------------

class BkgErrFiltParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(BkgErrFiltParameters, saber::SaberBlockParametersBase)
 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
  oops::RequiredParameter<float> oceanDepthMin{"ocean_depth_min", this};
  oops::RequiredParameter<double> rescaleBkgerr{"rescale_bkgerr", this};
  oops::RequiredParameter<double> efoldZ{"efold_z", this};
};

// ----------------------------------------------------------------------------------------

class BkgErrFilt : public saber::SaberOuterBlockBase {
 public:
  static const std::string classname() {return "soca::BkgErrFilt";}
  typedef BkgErrFiltParameters Parameters_;

  BkgErrFilt(const oops::GeometryData &,
    const oops::Variables &,
    const eckit::Configuration &,
    const Parameters_ &,
    const oops::FieldSet3D &,
    const oops::FieldSet3D &);

  virtual ~BkgErrFilt() = default;

  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override {}  // Not implemented

 private:
  void print(std::ostream &) const override;

  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;

  atlas::Field mult3D_;
  atlas::Field mult2D_;
  atlas::Field mask_;
};

// ----------------------------------------------------------------------------------------

}  // namespace soca
