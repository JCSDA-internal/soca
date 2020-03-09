/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_BALANCE_BALANCE_H_
#define SOCA_TRANSFORMS_BALANCE_BALANCE_H_

#include <ostream>
#include <string>

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

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

/// SOCA linear change of variable
class Balance: public util::Printable {
 public:
  struct Ftn{};

  static const std::string classname() {return "soca::Balance";}

  explicit Balance(const State &, const State &,
                   const Geometry &, const eckit::Configuration &);
  ~Balance();

/// Perform linear transforms
  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  void print(std::ostream &) const override;
  const State & traj_;
  Ftn * ftn_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_TRANSFORMS_BALANCE_BALANCE_H_
