
#ifndef SOCA_MODEL_MODELBIASINCREMENT_H_
#define SOCA_MODEL_MODELBIASINCREMENT_H_

#include <iostream>

#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace soca {
  class ModelBias;
  class ModelBiasCovariance;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasIncrement : public util::Printable {
 public:
/// Constructor, destructor
  ModelBiasIncrement(const Geometry &, const eckit::Configuration &) {}
  ModelBiasIncrement(const ModelBiasIncrement &, const bool) {}
  ModelBiasIncrement(const ModelBiasIncrement &, const eckit::Configuration &) {}
  ~ModelBiasIncrement() {}

/// Linear algebra operators
  void diff(const ModelBias &, const ModelBias &) {}
  void zero() {}
  ModelBiasIncrement & operator=(const ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator+=(const ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator-=(const ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator*=(const double) {return *this;}
  void axpy(const double, const ModelBiasIncrement &) {}
  double dot_product_with(const ModelBiasIncrement &) const {return 0.0;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  explicit ModelBiasIncrement(const ModelBiasCovariance &);
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_MODELBIASINCREMENT_H_
