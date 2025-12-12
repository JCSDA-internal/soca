/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/NumericConstraints.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

#include "soca/Utils/OceanSmoother.h"

namespace soca
{

/*
 * @brief A Saber Block to calculate the background error standard deviation of ocean fields
 * (temperature, unbalanced salinity, and unbalanced ssh) using a parametric
 * model based on the background vertical temperature gradient and mixed layer depth.
 *
 * @details For details see Weaver, A. T., Deltel, C., Machu, É., Ricci, S., & Daget, N.
 * (2005). A multivariate balance operator for variational ocean data
 * assimilation. Quarterly Journal of the Royal Meteorological Society: A journal
 * of the atmospheric sciences, applied meteorology and physical oceanography,
 * 131(613), 3605-3625.
 */
class ParametricOceanStdDev : public saber::SaberOuterBlockBase {
 public:
  // ----------------------------------------------------------------------------------------
  /// @brief Yaml parameters for ParametricOceanStdDev
  class Parameters : public saber::SaberBlockParametersBase {
    OOPS_CONCRETE_PARAMETERS(Parameters, saber::SaberBlockParametersBase)

   public:
    // --------------------------------------------------------------------------------------
    /// @brief yaml parameters that apply to all variables
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
    /// @brief yaml parameters for the temperature variable
    class Tocn : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(Tocn, BkgErrVariable)
     public:
      // special constructor to set default min and max
      Tocn() : BkgErrVariable("sea_water_potential_temperature", 0.1, 2.0) {}
      oops::Parameter<eckit::LocalConfiguration> sst{"sst",
        eckit::LocalConfiguration().set("fixed value", 1.0), this};
      oops::Parameter<double> dz{"dz", 20.0, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> efold{"efold", 500, this, {oops::minConstraint(0.0)}};
      oops::Parameter<double> minLayerThickness{"min layer thickness",
        "The minimum layer thickenss to use for calculating dt/dz, layers smaller than this"
        " are masked out", 0.01, this, {oops::minConstraint(0.0)}};
    };

    // --------------------------------------------------------------------------------------
    /// @brief yaml parameters for the unbalanced ssh variable
    class Ssh : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(Ssh, BkgErrVariable)
     public:
      // special constructor to set default min and max
      Ssh() : BkgErrVariable("sea_surface_height_above_geoid", 0.0, 0.1) {}
      oops::Parameter<double> phiEx{"phi ex", 20.0, this,
        {oops::minConstraint(0.0), oops::maxConstraint(90.0)}};
    };

    // --------------------------------------------------------------------------------------
    /// @brief yaml parameters for the unbalanced salinity variable
    class Socn : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS_ENABLE_COPY_AND_MOVE(Socn, BkgErrVariable)
     public:
      Socn() : BkgErrVariable("sea_water_salinity", 0.01, 0.25) {}  // set default min and max
      oops::Parameter<std::string> mldVariableName{"mld variable",
        "The name of the mixed layer depth variable in the background passed to the constructor",
        "ocean_mixed_layer_thickness", this};
      oops::Parameter<double> mldMax{"mld max", 400.0, this, {oops::minConstraint(0.0)}};
    };

    // --------------------------------------------------------------------------------------
    /// @brief yaml parameters for other variables that are not temperature, ssh, or salinity
    class OtherVar : public BkgErrVariable {
      OOPS_CONCRETE_PARAMETERS(OtherVar, BkgErrVariable)
     public:
      oops::RequiredParameter<std::string> varName{"variable", this};
      oops::Parameter<double> fractionOfBkg{"fraction of background", 0.0, this,
        {oops::minConstraint(0.0)}};
    };

    // --------------------------------------------------------------------------------------

    oops::Variables mandatoryActiveVars() const override { return oops::Variables(); }
    oops::Parameter<Tocn> tocn{"temperature", Tocn(), this};
    oops::Parameter<Socn> socn{"unbalanced salinity", Socn(), this};
    oops::Parameter<Ssh> ssh{"unbalanced ssh", Ssh(), this};
    oops::OptionalParameter<std::vector<OtherVar>> otherVars{"other variables", this};

    oops::OptionalParameter<OceanSmoother::Parameters> smoother{"smoother", this};
    oops::OptionalParameter<eckit::LocalConfiguration> saveDiags{"save diagnostics",
      "If present, save the calculated stddev to the given file",
      this};

    // other parameters that have defaults appropriate for the current variable
    // names in soca, but are here to keep me from hardcoding variable names.
    oops::Parameter<std::string> maskVariable{"mask variable",
      "The name of the geometry variable to use as a 2D horizontal land mask",
      "interp_mask", this};
    oops::Parameter<std::string> depthVariable{"depth variable",
      "The name of the depth state variable in the background passed to the constructor",
      "sea_water_depth", this};
  };

  // ----------------------------------------------------------------------------------------

  /// @brief This wrapper class is a bit of a hack. It allows us to do the yaml
  /// validation localy which otherwise is not catching extra erroneous parameters.
  class ParametersWrapper : public Parameters {
    OOPS_CONCRETE_PARAMETERS(ParametersWrapper, Parameters)
   public:
    oops::ConfigurationParameter fullConfig{this};
  };

  // ----------------------------------------------------------------------------------------

  typedef ParametersWrapper Parameters_;

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

}  // namespace soca
