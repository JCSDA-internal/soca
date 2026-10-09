/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

#include "soca/LinearVariableChange/Base/LinearVariableChangeBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace soca {
  class State;
  class Increment;
  class Geometry;
}

// -----------------------------------------------------------------------------

namespace soca {

/// SOCA balance: inverse
class InverseBalance: public LinearVariableChangeBase {
 public:
  static const std::string classname() {return "soca::Balance";}

  explicit InverseBalance(const State &, const State &,
                          const Geometry &, const eckit::Configuration &);
  ~InverseBalance();

/// Perform linear transforms
  void multiply(const Increment &, Increment &) const override;
  void multiplyAD(const Increment &, Increment &) const override;

 private:
  void print(std::ostream &) const override;
  int keyFtnConfig_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
