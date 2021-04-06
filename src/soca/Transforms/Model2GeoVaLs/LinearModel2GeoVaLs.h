/*
 * (C) Copyright 2021-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRANSFORMS_MODEL2GEOVALS_LINEARMODEL2GEOVALS_H_
#define SOCA_TRANSFORMS_MODEL2GEOVALS_LINEARMODEL2GEOVALS_H_

#include <string>
#include <memory>

#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
  class Geometry;
  class Increment;


class LinearModel2GeoVaLs: public util::Printable,
                           private util::ObjectCounter<LinearModel2GeoVaLs> {
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
  void print(std::ostream &) const override {}
};

}  // namespace soca

#endif  // SOCA_TRANSFORMS_MODEL2GEOVALS_LINEARMODEL2GEOVALS_H_
