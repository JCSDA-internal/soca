/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_BKGERR_BKGERR_H_
#define SOCA_TRANSFORMS_BKGERR_BKGERR_H_

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
class BkgErr: public util::Printable {
 public:
  static const std::string classname() {return "soca::BkgErr";}

  explicit BkgErr(const State &, const State &, const Geometry &,
                  const eckit::Configuration &);
  ~BkgErr();

  /// Read variance of background error, conf needs to be consistent with geom
  State initBkgErr(Geometry geom, const eckit::Configuration & conf) const
  {
    const eckit::LocalConfiguration lconf(conf, "background error variance");
    State bkgerr_variance(geom, lconf);
    return bkgerr_variance;
  }

/// Perform linear transforms
  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  void print(std::ostream &) const override;
  int keyFtnConfig_;
  State bkgerr_variance_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_TRANSFORMS_BKGERR_BKGERR_H_
