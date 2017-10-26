
#ifndef _MODEL_OBSFRACTION_H_
#define _MODEL_OBSFRACTION_H_

#include <ostream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "model/ObsSpace.h"
#include "model/Observation.h"
#include "model/ObsFractionTLAD.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace mom5cice5 {
  class Gom;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVec;

  // -----------------------------------------------------------------------------
  /// Sea-ice Fraction observation for model.
  /*!
   *  ObsFraction for model inherits from ObsEquivalent.
   */

  class ObsFraction : public Observation,
    private util::ObjectCounter<ObsFraction> {
  public:
      static const std::string classname() {return "mom5cice5::ObsFraction";}

      ObsFraction(ObsSpace &, const eckit::Configuration &);
      virtual ~ObsFraction();

      // Obs Operator
      void obsEquiv(const Gom &, ObsVec &, const ObsBias &) const;

      // Is there a way to put this in the TLAD class?
      LinearObsOp * getTLAD() const {return new ObsFractionTLAD(obsdb_, keyOperStrm_);}

      // Other
      void generateObsError(const eckit::Configuration &);
      boost::shared_ptr<const Variables> variables() const {return varin_;}

      int & toFortran() {return keyOperStrm_;}
      const int & toFortran() const {return keyOperStrm_;}

  private:
      void print(std::ostream &) const;
      ObsSpace & obsdb_;
      const std::string obsname_;
      int keyOperStrm_;
      boost::shared_ptr<const Variables> varin_;
    };
  // -----------------------------------------------------------------------------

}  // namespace mom5cice5
#endif  // _MODEL_OBSFRACTION_H_
