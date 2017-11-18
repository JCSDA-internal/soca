
#ifndef SOCA_MODEL_SOCAMODEL_H_
#define SOCA_MODEL_SOCAMODEL_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "model/Fortran.h"
#include "model/Geometry/Geometry.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
  //class F90traj;
  class ModelBias;
  class Fields;
  class State;

// -----------------------------------------------------------------------------
/// SOCA model definition.
/*!
 *  SOCA nonlinear model definition and configuration parameters.
 */

class Model: public util::Printable,
               private boost::noncopyable,
               private util::ObjectCounter<Model> {
 public:
  static const std::string classname() {return "soca::Model";}

  Model(const Geometry &, const eckit::Configuration &);
  ~Model();

/// Prepare model integration
  void initialize(State &) const;

/// Model integration
  void step(State &, const ModelBias &) const;
  int saveTrajectory(State &, const ModelBias &) const;

/// Finish model integration
  void finalize(State &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}

 private:
  void print(std::ostream &) const;
  int keyConfig_;
  util::Duration tstep_;
  const Geometry geom_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_MODEL_SOCAMODEL_H_
