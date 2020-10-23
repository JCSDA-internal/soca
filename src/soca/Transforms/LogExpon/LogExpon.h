/*
 * (C) Copyright 2020-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_LOGEXPON_LOGEXPON_H_
#define SOCA_TRANSFORMS_LOGEXPON_LOGEXPON_H_

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

class LogExpon: public util::Printable {
 public:
  static const std::string classname() {return "soca::LogExpon";}

  explicit LogExpon(const Geometry &,
                     const eckit::Configuration &);
  ~LogExpon();

  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

  std::vector<std::string> initTrans(const eckit::Configuration & conf,
                                      const std::string & trvar) const
  {
    return conf.getStringVector("Transform."+trvar);
  }

 private:
  void print(std::ostream &) const override;
  const oops::Variables trvars_;
};
// -----------------------------------------------------------------------------
}  // namespace soca

#endif  // SOCA_TRANSFORMS_LOGEXPON_LOGEXPON_H_
