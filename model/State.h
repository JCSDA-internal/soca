
#ifndef MOM5CICE5_MODEL_MOM5CICE5STATE_H_
#define MOM5CICE5_MODEL_MOM5CICE5STATE_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "model/Fields.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mom5cice5 {
  class Geometry;
  class Increment;
  class Variables;

/// MOM5CICE5 model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class State : public util::Printable,
                private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "mom5cice5::State";}

/// Constructor, destructor
  State(const Geometry &, const Variables &, const util::DateTime &);  // Is it used?
  State(const Geometry &, const eckit::Configuration &);
  State(const Geometry &, const State &);
  State(const State &);
  virtual ~State();
  State & operator=(const State &);

/// Interpolate to observation location
//  void interpolate(const LocQG &, GomQG &) const;

/// Interpolate full fields
///  void changeResolution(const State & xx);

/// Interactions with Increment
  State & operator+=(const Increment &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}

/// Access to fields
  Fields & fields() {return *fields_;}
  const Fields & fields() const {return *fields_;}

  boost::shared_ptr<const Geometry> geometry() const {
    return fields_->geometry();
  }

/// Other
  void activateModel();
  void deactivateModel();

  void zero();
  void accumul(const double &, const State &);

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<Fields> fields_;
  boost::scoped_ptr<Fields> stash_;
};
// -----------------------------------------------------------------------------

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_MOM5CICE5STATE_H_
