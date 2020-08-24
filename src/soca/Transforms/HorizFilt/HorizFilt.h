/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_HORIZFILT_HORIZFILT_H_
#define SOCA_TRANSFORMS_HORIZFILT_HORIZFILT_H_

#include <ostream>
#include <memory>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace soca {
  class Fields;
  class State;
  class Increment;
  class Geometry;
}

// -----------------------------------------------------------------------------

namespace soca {

/// SOCA linear change of variable
class HorizFilt: public util::Printable {
 public:
  static const std::string classname() {return "soca::HorizFilt";}

  explicit HorizFilt(const State &, const State &, const Geometry &,
                  const eckit::Configuration &);
  ~HorizFilt();

/// Perform linear transforms
  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  void print(std::ostream &) const override;
  int keyFtnConfig_;
  std::unique_ptr<const Geometry> geom_;
  oops::Variables vars_;
  const State traj_;
  unsigned int niter_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_TRANSFORMS_HORIZFILT_HORIZFILT_H_
