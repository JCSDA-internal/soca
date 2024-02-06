/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>
#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace soca {
  class Geometry;
}

// -----------------------------------------------------------------------------

namespace soca {

class ModelData : public util::Printable {
 public:
  static const std::string classname() {return "soca::ModelData";}
  static const oops::Variables defaultVariables() {
    return oops::Variables(std::vector<std::string>({"surface_temperature_where_sea"}));
  }

  explicit ModelData(const Geometry &) {}
  ~ModelData() {}

  const eckit::LocalConfiguration modelData() const {return eckit::LocalConfiguration();}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace soca
