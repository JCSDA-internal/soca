
#ifndef MOM5CICE5_MODEL_LINEAROBSOP_H_
#define MOM5CICE5_MODEL_LINEAROBSOP_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "model/ObsOperator/Observation.h"
#include "eckit/config/Configuration.h"

namespace util {
  class DateTime;
}

namespace mom5cice5 {
  class Gom;
  class ObsSpace;
  class ObsVec;
  class ObsBias;
  class ObsBiasIncrement;
  class Variables;

  // -----------------------------------------------------------------------------

  class LinearObsOp : private boost::noncopyable {
  public:
    static LinearObsOp * create(const Observation & obs) {return obs.getTLAD();}

    LinearObsOp() {}
    virtual ~LinearObsOp() {}

    // Obs Operators
    virtual void setTrajectory(const Gom &, const ObsBias &) =0;
    virtual void obsEquivTL(const Gom &, ObsVec &, const ObsBiasIncrement &) const =0;
    virtual void obsEquivAD(Gom &, const ObsVec &, ObsBiasIncrement &) const =0;

    // Other
    virtual boost::shared_ptr<const Variables> variables() const =0;

    virtual int & toFortran() =0;
    virtual const int & toFortran() const =0;
  };

  // -----------------------------------------------------------------------------

}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_LINEAROBSOP_H_
