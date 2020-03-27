/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_ANA2MODEL_ANA2MODEL_H_
#define SOCA_TRANSFORMS_ANA2MODEL_ANA2MODEL_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "soca/Geometry/Geometry.h"
#include "oops/util/Printable.h"
#include "oops/base/Variables.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// SOCA nonlinear change of variable

class Ana2Model: public util::Printable {
 public:
  static const std::string classname() {return "soca::Ana2Model";}

  explicit Ana2Model(const Geometry &,
                     const eckit::Configuration &);
  ~Ana2Model();

  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

  std::vector<std::string> initRotate(const eckit::Configuration & conf,
                                      const std::string & uv) const
  {
    return conf.getStringVector("rotate."+uv);
  }

 private:
  void print(std::ostream &) const override;
  const oops::Variables uvars_;
  const oops::Variables vvars_;
};
// -----------------------------------------------------------------------------
}  // namespace soca

#endif  // SOCA_TRANSFORMS_ANA2MODEL_ANA2MODEL_H_
