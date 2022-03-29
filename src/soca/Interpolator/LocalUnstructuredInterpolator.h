/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/generic/InterpolatorUnstructured.h"
#include "oops/util/Printable.h"


namespace eckit {
  class Configuration;
  class LocalConfiguration;
}

namespace oops {
  class Variables;
}

namespace soca {
  class Geometry;
  class Increment;
  class State;
}


namespace soca {

// -----------------------------------------------------------------------------

class LocalUnstructuredInterpolator : public util::Printable {
 public:
  LocalUnstructuredInterpolator(const eckit::Configuration &, const Geometry &,
                                const std::vector<double> &);
  ~LocalUnstructuredInterpolator() {}

  void apply(const oops::Variables &, const State &, std::vector<double> &) const;
  void apply(const oops::Variables &, const Increment &, std::vector<double> &) const;
  void applyAD(const oops::Variables &, Increment &, const std::vector<double> &) const;

 private:
  const std::shared_ptr<oops::InterpolatorUnstructured> getInterpolator(const std::string &) const;
  void print(std::ostream &) const;

  mutable std::shared_ptr<oops::InterpolatorUnstructured> interp_[6];

  const std::shared_ptr<const Geometry> geom_;
  const eckit::LocalConfiguration config_;
  const std::vector<double> latlon_out_;
};

}  // namespace soca
