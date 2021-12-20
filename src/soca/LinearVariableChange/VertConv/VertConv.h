/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "soca/State/State.h"

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

#include "soca/LinearVariableChange/Base/LinearVariableChangeBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace soca {
  class Geometry;
  class Increment;
}

// -----------------------------------------------------------------------------

namespace soca {

/// SOCA linear change of variable
class VertConv: public LinearVariableChangeBase {
 public:
  static const std::string classname() {return "soca::VertConv";}

  explicit VertConv(const State &, const State &, const Geometry &,
                    const eckit::Configuration &);
  ~VertConv();

/// Perform linear transforms
  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  void print(std::ostream &) const override;
  int keyFtnConfig_;
  const State bkg_lr_;
  const Geometry geom_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
