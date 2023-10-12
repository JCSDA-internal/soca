/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"


// --------------------------------------------------------------------------------------
// Forward declarations
namespace soca {
  class Geometry;
}

// --------------------------------------------------------------------------------------

namespace soca {

// --------------------------------------------------------------------------------------

class ExplicitDiffusionIOParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionIOParameters, oops::Parameters)
 public:
  oops::RequiredParameter<std::string> filename{"filename", this};
};

// --------------------------------------------------------------------------------------

class ExplicitDiffusionCalibrationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionCalibrationParameters, oops::Parameters)
 public:
  // TODO, formalize
  oops::RequiredParameter<eckit::LocalConfiguration> normalization{"normalization", this};
  oops::RequiredParameter<eckit::LocalConfiguration> scales{"scales", this};
  oops::OptionalParameter<ExplicitDiffusionIOParameters> write {"write", this};
};

// --------------------------------------------------------------------------------------

class ExplicitDiffusionParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ExplicitDiffusionParameters, saber::SaberBlockParametersBase)
 public:
  oops::OptionalParameter<ExplicitDiffusionIOParameters> read{"read", this};
  oops::OptionalParameter<ExplicitDiffusionCalibrationParameters> 
    calibration{"calibration", this};  
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};

  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// --------------------------------------------------------------------------------------

class ExplicitDiffusion : public saber::SaberCentralBlockBase {
 public:
  static const std::string classname() { return "soca::ExplicitDiffusion"; }
  typedef ExplicitDiffusionParameters Parameters_;

  ExplicitDiffusion(const oops::GeometryData &,
                    const oops::Variables &,
                    const eckit::Configuration &,
                    const Parameters_ &,
                    const oops::FieldSet3D &,
                    const oops::FieldSet3D &);

  void randomize(atlas::FieldSet &) const override;
  void multiply(atlas::FieldSet &) const override;

  void directCalibration(const std::vector<atlas::FieldSet> &) override;
  void read() override;  
  void write() const override;

 private:
  void print(std::ostream &) const override;
  std::shared_ptr<Geometry> geom_;
  int keyFortran_;
  oops::Variables vars_;
  Parameters_ params_;
};

// --------------------------------------------------------------------------------------

}  // namespace soca