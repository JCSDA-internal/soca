/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/FieldSets.h"
#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/NumericConstraints.h"

#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"

// -----------------------------------------------------------------------------------------
namespace oops {
  class FieldSet3D;
}
namespace eckit {
  class Configuration;
}

// -----------------------------------------------------------------------------------------

namespace soca
{

// -----------------------------------------------------------------------------------------
// Parameters
class ParametricOceanStdDevBound : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ParametricOceanStdDevBound, oops::Parameters)
 public:
  oops::Parameter<double> min{"min", 0.0, this, {oops::minConstraint(0.0)} };
  oops::Parameter<double> max{"max", std::numeric_limits<double>::max(), this, {oops::minConstraint(0.0)} };
};

class ParametricOceanStdDevTocn : public ParametricOceanStdDevBound {
  OOPS_CONCRETE_PARAMETERS(ParametricOceanStdDevTocn, ParametricOceanStdDevBound)
 public:
  oops::RequiredParameter<std::string> sstFile{"sst file", this};
  oops::RequiredParameter<double> dz{"dz", this, {oops::minConstraint(0.0)}};
  oops::RequiredParameter<double> efold{"efold", this, {oops::minConstraint(0.0)}};
};

class ParametricOceanStdDevSsh : public ParametricOceanStdDevBound {
  OOPS_CONCRETE_PARAMETERS(ParametricOceanStdDevSsh, ParametricOceanStdDevBound)
 public:
  oops::RequiredParameter<double> phi_ex{"phi_ex", this, 
    {oops::minConstraint(0.0), oops::maxConstraint(90.0)}};
};

class ParametricOceanStdDevOther : public ParametricOceanStdDevBound {
  OOPS_CONCRETE_PARAMETERS(ParametricOceanStdDevOther, ParametricOceanStdDevBound)
 public:
  oops::RequiredParameter<double> fractionOfBkg{"fraction of bkg", this, 
    {oops::minConstraint(0.0)}};
};

class ParametricOceanStdDevParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(ParametricOceanStdDevParameters, saber::SaberBlockParametersBase)
 public:
  oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
  oops::RequiredParameter<ParametricOceanStdDevTocn> tocn{"temperature", this};
  oops::RequiredParameter<ParametricOceanStdDevBound> socn{"unbalanced salinity", this};
  oops::RequiredParameter<ParametricOceanStdDevSsh> ssh{"unbalanced ssh", this};
  oops::OptionalParameter<std::map<std::string, ParametricOceanStdDevOther> > others{"others", this};
  oops::OptionalParameter<eckit::LocalConfiguration> output{"output", this};
  oops::OptionalParameter<eckit::LocalConfiguration> smoother{"smoother", this};
};

// -----------------------------------------------------------------------------------------

class ParametricOceanStdDev : public saber::SaberOuterBlockBase {
 public:
  static const std::string classname() {return "soca::ParametricOceanStdDev";}

  typedef ParametricOceanStdDevParameters Parameters_;

  ParametricOceanStdDev(const oops::GeometryData &,
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
  oops::FieldSet3D bkgErr_;
};

}  // namespace soca
