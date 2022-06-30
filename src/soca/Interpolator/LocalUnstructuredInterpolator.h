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

#include "oops/util/Printable.h"

namespace atlas {
  class FieldSet;
}

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
  class UnstructuredInterpolator;
}


namespace soca {

// -----------------------------------------------------------------------------

class LocalUnstructuredInterpolator : public util::Printable {
 public:
  LocalUnstructuredInterpolator(const eckit::Configuration &, const Geometry &,
                                const std::vector<double> &, const std::vector<double> &);
  ~LocalUnstructuredInterpolator() {}

  void apply(const oops::Variables &, const atlas::FieldSet &, const std::vector<bool> &,
             std::vector<double> &) const;
  void applyAD(const oops::Variables &, atlas::FieldSet &, const std::vector<bool> &,
               const std::vector<double> &) const;

 private:
  const std::shared_ptr<UnstructuredInterpolator> getInterpolator(const std::string &) const;
  void print(std::ostream &) const;

  std::shared_ptr<UnstructuredInterpolator> interp_[6];
  const Geometry & geom_;
};

}  // namespace soca
