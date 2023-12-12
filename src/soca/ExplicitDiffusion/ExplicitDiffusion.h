/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "soca/ExplicitDiffusion/ExplicitDiffusionParameters.h"
#include "saber/blocks/SaberCentralBlockBase.h"


// --------------------------------------------------------------------------------------
// Forward declarations
namespace soca {
  class Geometry;
}

// --------------------------------------------------------------------------------------

namespace soca {

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

  void randomize(oops::FieldSet3D &) const override;
  void multiply(oops::FieldSet3D &) const override;

  void directCalibration(const std::vector<oops::FieldSet3D> &) override;
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
