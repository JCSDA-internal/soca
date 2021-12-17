/*
 * (C) Copyright 2021-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/util/Printable.h"

#include "soca/LinearVariableChange/Base/LinearVariableChangeBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
  class Geometry;
  class Increment;


class LinearModel2GeoVaLs: public LinearVariableChangeBase {
 public:
  static const std::string classname() {return "soca::LinearModel2GeoVaLs";}

  explicit LinearModel2GeoVaLs(const State &, const State &, const Geometry &,
                                  const eckit::Configuration &);
  ~LinearModel2GeoVaLs();

  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  std::unique_ptr<const Geometry> geom_;
  void print(std::ostream &) const;
};

}  // namespace soca
