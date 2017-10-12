
#ifndef MOM5CICE5_MODEL_OBSERVATION_H_
#define MOM5CICE5_MODEL_OBSERVATION_H_

#include <map>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "util/Printable.h"

namespace util {
  class DateTime;
}

namespace mom5cice5 {
  class Gom;
  class LinearObsOp;
  class ObsSpace;
  class ObsVec;
  class ObsBias;
  class Variables;

/// Observations for MOM5CICE5 model.
/*!
 *  Base class for MOM5CICE5 observations.
 */

// -----------------------------------------------------------------------------

class Observation;

/// ObsFactory
class ObsFactory {
 public:
  static Observation * create(ObsSpace &, const eckit::Configuration &);
  virtual ~ObsFactory() { makers_->clear(); }
 protected:
  explicit ObsFactory(const std::string &);
 private:
  virtual Observation * make(ObsSpace &, const eckit::Configuration &) =0;
  static std::map < std::string, ObsFactory * > * makers_;
};

template<class T>
class ObsMaker : public ObsFactory {
  virtual Observation * make(ObsSpace & odb, const eckit::Configuration & c) {return new T(odb, c);}
 public:
  explicit ObsMaker(const std::string & name) : ObsFactory(name) {}
};

// -----------------------------------------------------------------------------

class Observation : public util::Printable,
                      private boost::noncopyable {
 public:
  static Observation * create(ObsSpace & odb, const eckit::Configuration & conf)
    {return ObsFactory::create(odb, conf);}

  virtual ~Observation() {}

// Obs Operators
  virtual void obsEquiv(const Gom &, ObsVec &, const ObsBias &) const =0;

// Get TLAD obs operator (should be more like static create(...)?
  virtual LinearObsOp * getTLAD() const =0;

// Other
  virtual void generateObsError(const eckit::Configuration &) =0;
  virtual boost::shared_ptr<const Variables> variables() const =0;

  virtual int & toFortran() =0;
  virtual const int & toFortran() const =0;

 private:
  virtual void print(std::ostream &) const =0;
};

// -----------------------------------------------------------------------------

}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_OBSERVATION_H_
