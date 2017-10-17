
#ifndef MOM5CICE5_MODEL_MODELBIAS_H_
#define MOM5CICE5_MODEL_MODELBIAS_H_

#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mom5cice5 {
  class Geometry;
  class ModelBiasIncrement;

/// Model error for the MOM5CICE5 model.
/*!
 * This class is used to manipulate parameters of the model that
 * can be estimated in the assimilation. This includes model bias for
 * example but could be used for other parameters to be estimated.
 * This is sometimes referred to as augmented state or augmented
 * control variable in the litterature.
 * The augmented state is understood here as an augmented 4D state.
 */

// -----------------------------------------------------------------------------

class ModelBias : public util::Printable,
                  private boost::noncopyable,
                  private util::ObjectCounter<ModelBias> {
 public:
  static const std::string classname() {return "mom5cice5::ModelBias";}

  ModelBias(const Geometry &, const eckit::Configuration &) {}
  ModelBias(const Geometry &, const ModelBias &) {}
  ModelBias(const ModelBias &, const bool) {}
  ~ModelBias() {}

  ModelBias & operator+=(const ModelBiasIncrement &) {return *this;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_MODELBIAS_H_
