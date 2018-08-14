/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_CHANGEVAR_H_
#define SOCA_MODEL_CHANGEVAR_H_

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
  class Increment;

// -----------------------------------------------------------------------------
/// SOCA linear change of variable

class Kst: public util::Printable {
 public:
  static const std::string classname() {return "soca::Kst";}

  explicit Kst(const State &, const State &, const eckit::Configuration &);
  ~Kst();

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
#endif  // SOCA_MODEL_CHANGEVAR_H_
