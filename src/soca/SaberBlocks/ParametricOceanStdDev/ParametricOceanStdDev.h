/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <limits>

#include "oops/util/parameters/NumericConstraints.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

namespace soca
{

class ParametricOceanStdDev : public saber::SaberOuterBlockBase {
 public:
  // ----------------------------------------------------------------------------------------
  class Parameters : public saber::SaberBlockParametersBase {
    OOPS_CONCRETE_PARAMETERS(Parameters, saber::SaberBlockParametersBase)
   public:
    // --------------------------------------------------------------------------------------
    class Bounds : public oops::Parameters {
      OOPS_CONCRETE_PARAMETERS(Bounds, oops::Parameters)
     public:
      oops::Parameter<double> min{"min", 0.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> max{"max", std::numeric_limits<double>::max(),
        this, {oops::minConstraint(0.0)}};
    };

    // --------------------------------------------------------------------------------------
    class Tocn : public Bounds {
      OOPS_CONCRETE_PARAMETERS(Tocn, Bounds)
     public:
    };

    // --------------------------------------------------------------------------------------
    class Socn : public Bounds {
      OOPS_CONCRETE_PARAMETERS(Socn, Bounds)
     public:
    };

    // --------------------------------------------------------------------------------------
    class Ssh : public Bounds {
      OOPS_CONCRETE_PARAMETERS(Ssh, Bounds)
     public:
    };

    // --------------------------------------------------------------------------------------

    oops::Variables mandatoryActiveVars() const override { return oops::Variables(); }
    oops::Parameter<Tocn> tocn{"temperature", Tocn(), this};
    oops::Parameter<Socn> socn{"unbalanced salinity", Socn(), this};
    oops::Parameter<Ssh> ssh{"unbalanced ssh", Ssh(), this};

  };
  // ----------------------------------------------------------------------------------------

  typedef Parameters Parameters_;

  explicit ParametricOceanStdDev(const oops::GeometryData &,
                                  const oops::Variables &,
                                  const eckit::Configuration &,
                                  const Parameters_ &,
                                  const oops::FieldSet3D &,
                                  const oops::FieldSet3D &);

  virtual ~ParametricOceanStdDev() = default;

  const oops::GeometryData & innerGeometryData() const override { return innerGeometryData_; }
  const oops::Variables & innerVars() const override { return innerVars_; }

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;

  const oops::GeometryData & innerGeometryData_;
  oops::Variables innerVars_;
};

} // namespace soca
