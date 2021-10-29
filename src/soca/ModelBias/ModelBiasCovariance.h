/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODELBIAS_MODELBIASCOVARIANCE_H_
#define SOCA_MODELBIAS_MODELBIASCOVARIANCE_H_

#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace soca {
  class Geometry;
  class ModelBias;
  class ModelBiasIncrement;
}

// -----------------------------------------------------------------------------

namespace soca {

class ModelBiasCovariance : public util::Printable,
                            private boost::noncopyable,
                            private util::ObjectCounter<ModelBiasCovariance> {
 public:
  static const std::string classname() {return "soca::ModelBiasCovariance";}

/// Constructor, destructor
  ModelBiasCovariance(const eckit::Configuration & conf,
                      const Geometry &): conf_(conf) {}
  ~ModelBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ModelBias &, const Geometry &) {}
  void multiply(const ModelBiasIncrement &, ModelBiasIncrement &) const {}
  void inverseMultiply(const ModelBiasIncrement &,
                       ModelBiasIncrement &) const {}
  void randomize(ModelBiasIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODELBIAS_MODELBIASCOVARIANCE_H_
