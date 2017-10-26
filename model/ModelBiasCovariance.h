
#ifndef MOM5CICE5_MODEL_MODELBIASCOVARIANCE_H_
#define MOM5CICE5_MODEL_MODELBIASCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace mom5cice5 {
  class ModelBias;
  class ModelBiasIncrement;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasCovariance : public util::Printable,
                            private boost::noncopyable,
                            private util::ObjectCounter<ModelBiasCovariance> {
 public:
  static const std::string classname() {return "mom5cice5::ModelBiasCovariance";}

/// Constructor, destructor
  ModelBiasCovariance(const eckit::Configuration & conf, const Geometry &): conf_(conf) {}
  ~ModelBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ModelBias &, const Geometry &) {}
  void multiply(const ModelBiasIncrement &, ModelBiasIncrement) const {}
  void inverseMultiply(const ModelBiasIncrement &, ModelBiasIncrement) const {}
  void randomize(ModelBiasIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_MODELBIASCOVARIANCE_H_
