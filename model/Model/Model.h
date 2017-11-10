
#ifndef MOM5CICE5_MODEL_MOM5CICE5MODEL_H_
#define MOM5CICE5_MODEL_MOM5CICE5MODEL_H_

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

namespace mom5cice5 {
  //class F90traj;
  class ModelBias;
  class Fields;
  class State;

// -----------------------------------------------------------------------------
/// MOM5CICE5 model definition.
/*!
 *  MOM5CICE5 nonlinear model definition and configuration parameters.
 */

class Model: public util::Printable,
               private boost::noncopyable,
               private util::ObjectCounter<Model> {
 public:
  static const std::string classname() {return "mom5cice5::Model";}

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

}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_MOM5CICE5MODEL_H_
