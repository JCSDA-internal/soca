/*
 * (C) Copyright 2021-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

#include "soca/VariableChange/Base/VariableChangeBase.h"

namespace soca {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------

class VariableChangeParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters,
                           oops::VariableChangeParametersBase)
 public:
  // Wrapper to VariableChange parameters
  VariableChangeParametersWrapper variableChangeParametersWrapper{this};
};

// -----------------------------------------------------------------------------

class VariableChange : public util::Printable {
 public:
  static const std::string classname() {return "soca::VariableChange";}

  typedef VariableChangeParameters Parameters_;

  explicit VariableChange(const Parameters_ &, const Geometry &);
  ~VariableChange();

  void changeVar(State &, const oops::Variables &) const;
  void changeVarInverse(State &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<VariableChangeBase> variableChange_;
};

// -----------------------------------------------------------------------------

}  // namespace soca
