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

#include "eckit/config/Configuration.h"
#include "soca/Geometry/Geometry.h"
#include "oops/util/Printable.h"

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

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------
}  // namespace soca

#endif  // SOCA_TRANSFORMS_ANA2MODEL_ANA2MODEL_H_
