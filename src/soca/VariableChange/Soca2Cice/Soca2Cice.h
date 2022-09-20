/*
 * (C) Copyright 2022-2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "soca/Traits.h"

#include "soca/VariableChange/Base/VariableChangeBase.h"

namespace soca {

class Soca2Cice: public VariableChangeBase {
 public:
  const std::string classname() {return "soca::Soca2Cice";}

  Soca2Cice(const Geometry &, const eckit::Configuration &);
  ~Soca2Cice();

  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  const Geometry & geom_;
  int keySoca2Cice_;
  void print(std::ostream &) const override {}
};

}  // namespace soca
