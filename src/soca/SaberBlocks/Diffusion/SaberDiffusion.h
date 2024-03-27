/*
 * (C) Copyright 2023-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#pragma once

#include "saber/blocks/SaberCentralBlockBase.h"

#include "soca/SaberBlocks/Diffusion/SaberDiffusionParameters.h"

//forward declarations
namespace soca {
  class Diffusion;
}


namespace soca {

class SaberDiffusion : public saber::SaberCentralBlockBase {
 public:
  static const std::string classname() { return "soca::SaberDiffusion"; }
  typedef SaberDiffusionParameters Parameters_;

  SaberDiffusion(const oops::GeometryData &,
                 const oops::Variables &,
                 const eckit::Configuration &,
                 const Parameters_ &,
                 const oops::FieldSet3D &,
                 const oops::FieldSet3D &);

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;
  
  void read() override;
  void directCalibration(const oops::FieldSets &) override;

 private:
  void print(std::ostream &) const override {}
  Parameters_ params_;
  const oops::GeometryData & geom_;
  std::unique_ptr<Diffusion> diffusion_;

};

}