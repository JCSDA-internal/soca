/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_OBSBIASCOVARIANCE_H_
#define SOCA_MODEL_OBSBIASCOVARIANCE_H_

#include <ostream>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace soca {
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------

class ObsBiasCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsBiasCovariance> {
 public:
  static const std::string classname() {return "soca::ObsBiasCovariance";}

/// Constructor, destructor
  explicit ObsBiasCovariance(const eckit::Configuration &);
  ~ObsBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ObsBias &) {}
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void inverseMultiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;
  void randomize(ObsBiasIncrement &) const;

  const eckit::Configuration & config() const {return conf_;}
  bool active(const unsigned int ii) const {return variance_[ii] > 0.0;}

 private:
  void print(std::ostream &) const;
  const eckit::LocalConfiguration conf_;
  std::vector<double> variance_;
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_OBSBIASCOVARIANCE_H_
