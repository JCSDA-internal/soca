/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Traits.h"

#include "soca/VariableChange/Base/VariableChangeBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// SOCA nonlinear change of variable

class Ana2Model: public VariableChangeBase {
 public:
  const std::string classname() {return "soca::Ana2Model";}

  Ana2Model(const Geometry &, const eckit::Configuration &);
  ~Ana2Model();

  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

  std::vector<std::string> initRotate(const eckit::Configuration & conf,
                                      const std::string & uv) const
  {
    return conf.getStringVector("rotate."+uv);
  }

  std::vector<std::string> initTrans(const eckit::Configuration & conf,
                                      const std::string & trvar) const
  {
    return conf.getStringVector("log."+trvar);
  }

 private:
  void print(std::ostream &) const override;
  const oops::Variables uvars_;
  const oops::Variables vvars_;
  const oops::Variables logvars_;
};
// -----------------------------------------------------------------------------
}  // namespace soca
