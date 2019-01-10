/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_TRANSFORMS_VERTCONV_VERTCONV_H_
#define SRC_TRANSFORMS_VERTCONV_VERTCONV_H_

#include <ostream>
#include <string>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"
#include "eckit/config/Configuration.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
  class State;
  class Geometry;
  class Increment;

// -----------------------------------------------------------------------------
/// SOCA linear change of variable

class VertConv: public util::Printable {
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
  const State & traj_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SRC_TRANSFORMS_VERTCONV_VERTCONV_H_
