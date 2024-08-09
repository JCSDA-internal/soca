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

#include "soca/Utils/OceanSmoother.h"

namespace soca
{

class ParametricOceanStdDev : public saber::SaberOuterBlockBase {
 public:
  // ----------------------------------------------------------------------------------------
  // Yaml parameters for ParametricOceanStdDev
  class Parameters : public saber::SaberBlockParametersBase {
    OOPS_CONCRETE_PARAMETERS(Parameters, saber::SaberBlockParametersBase)
   public:
    // --------------------------------------------------------------------------------------
    class BkgErrVariable : public oops::Parameters {
      OOPS_ABSTRACT_PARAMETERS(BkgErrVariable, oops::Parameters)
     public:
      // a special constructor to set the default min and max
      BkgErrVariable(std::string defaultVarName, double defaultMin, double defaultMax)
        : oops::Parameters(),
          varName("variable", defaultVarName, this),
          min("min", defaultMin, this, {oops::minConstraint(0.0)}),
          max("max", defaultMax, this, {oops::minConstraint(0.0)}) {}
      oops::Parameter<std::string> varName{"variable", "", this};
      oops::Parameter<double> min{"min", 0.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> max{"max", std::numeric_limits<double>::max(),
        this, {oops::minConstraint(0.0)}};
      oops::Parameter<bool> smooth{"smooth", true, this};
    };

    // --------------------------------------------------------------------------------------
    class Tocn : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(Tocn, BkgErrVariable)
     public:
      Tocn() : BkgErrVariable("tocn", 0.1, 2.0) {}  // set default min and max
      oops::Parameter<eckit::LocalConfiguration> sst{"sst",
        eckit::LocalConfiguration().set("fixed value", 1.0), this};
      oops::Parameter<double> dz{"dz", 20.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> efold{"efold", 500, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> minLayerThickness{"min layer thickness",
        "The minimum layer thickenss to use for calculating dt/dz, layers smaller than this"
        " are masked out", 0.01, this, {oops::minConstraint(0.0)}};
    };

    // --------------------------------------------------------------------------------------
    class Ssh : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(Ssh, BkgErrVariable)
     public:
      Ssh() : BkgErrVariable("ssh", 0.0, 0.1) {}  // set default min and max
      oops::Parameter<double> phiEx{"phi ex", 20.0, this,
        {oops::minConstraint(0.0), oops::maxConstraint(90.0)}};
    };

    // --------------------------------------------------------------------------------------
    class Socn : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(Socn, BkgErrVariable)
     public:
      Socn() : BkgErrVariable("socn", 0.0, 0.25) {}  // set default min and max
      oops::Parameter<std::string> mldVariableName{"mld variable",
        "The name of the mixed layer depth variable in the background passed to the constructor",
        "mld", this};
      oops::Parameter<double> mldMax{"mld max", 400.0, this, {oops::minConstraint(0.0)}};
    };

    // --------------------------------------------------------------------------------------
    class OtherVar : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS(OtherVar, BkgErrVariable)
     public:
      oops::Parameter<double> fractionOfBkg{"fraction of background", 0.0, this,
        {oops::minConstraint(0.0)}};
    };

    // --------------------------------------------------------------------------------------

    oops::Variables mandatoryActiveVars() const override { return oops::Variables(); }
    oops::Parameter<Tocn> tocn{"temperature", Tocn(), this};
    oops::Parameter<Socn> socn{"unbalanced salinity", Socn(), this};
    oops::Parameter<Ssh> ssh{"unbalanced ssh", Ssh(), this};
    oops::OptionalParameter<std::map<std::string, OtherVar>> otherVars{"other variables", this};

    oops::OptionalParameter<OceanSmoother::Parameters> smoother{"smoother", this};
    oops::OptionalParameter<eckit::LocalConfiguration> saveDiags{"save diagnostics",
      "If present, save the calculated stddev to the given file",
      this};

    // other parameters that have defaults appropriate for the current variable
    // names in soca, but are here to keep me from hardcoding variable names.
    oops::Parameter<std::string> maskVariable{"mask variable",
      "The name of the geometry variable to use as a 2D horizontal land mask",
      "interp_mask", this};
    oops::Parameter<std::string> thicknessVariable{"thickness variable",
      "The name of the layer thickness state variable in the background passed to the constructor",
      "hocn", this};

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

  const oops::GeometryData & innerGeometryData() const override { return geom_; }
  const oops::Variables & innerVars() const override { return innerVars_; }

  void multiply(oops::FieldSet3D &) const override;
  void multiplyAD(oops::FieldSet3D &) const override;
  void leftInverseMultiply(oops::FieldSet3D &) const override;

 private:
  void print(std::ostream &) const override;
  void commonMultiply(oops::FieldSet3D &) const;

  const oops::GeometryData & geom_;
  oops::Variables innerVars_;
  atlas::FieldSet bkgErr_;
};

} // namespace soca
