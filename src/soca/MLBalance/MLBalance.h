/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "soca/MLBalance/MLBalanceParameters.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"


// --------------------------------------------------------------------------------------
// Forward declarations
namespace soca {
  class Geometry;
}

// --------------------------------------------------------------------------------------

namespace soca {

// --------------------------------------------------------------------------------------

class MLBalance : public saber::SaberOuterBlockBase {
 public:
  static const std::string classname() { return "soca::MLBalance"; }
  typedef MLBalanceParameters Parameters_;

  MLBalance(const oops::GeometryData &,
            const oops::Variables &,
            const eckit::Configuration &,
            const Parameters_ &,
            const oops::FieldSet3D &,
            const oops::FieldSet3D &);
  virtual ~MLBalance();

  const oops::Variables & innerVars() const override {return innerVars_;}
  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;

  void setupJac(const oops::FieldSet3D &,
                const eckit::mpi::Comm &,
                const eckit::Configuration &);

 private:
  void print(std::ostream &) const override;
  oops::Variables activeVars_;
  oops::Variables innerVars_;
  const oops::GeometryData & innerGeometryData_;
  Parameters_ params_;
  atlas::FieldSet jac_;
};

// --------------------------------------------------------------------------------------

}  // namespace soca
