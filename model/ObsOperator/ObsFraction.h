
#ifndef _MODEL_OBSFRACTION_H_
#define _MODEL_OBSFRACTION_H_

#include <ostream>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/interface/ObsOperatorBase.h"
#include "model/ObsSpace/ObsSpace.h"
//#include "model/ObsOperator/Observation.h"
//#include "model/ObsOperator/ObsFractionTLAD.h"
#include "util/ObjectCounter.h"
#include "model/Traits.h"

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

  class ObsFraction : public oops::ObsOperatorBase<Traits>,
                      private util::ObjectCounter<ObsFraction> {
  public:
      static const std::string classname() {return "mom5cice5::ObsFraction";}

      ObsFraction(const ObsSpace &, const eckit::Configuration &);
      virtual ~ObsFraction();

      // Obs Operator
      void obsEquiv(const Gom &, ObsVec &, const ObsBias &) const;

      // Other
      boost::shared_ptr<const Variables> variables() const {return varin_;}

      int & toFortran() {return keyOperFraction_;}
      const int & toFortran() const {return keyOperFraction_;}

  private:
      void print(std::ostream &) const;
      //const ObsSpace & obsdb_;
      const std::string obsname_;
      //int keyOperFraction_;
      F90hop keyOperFraction_;
      boost::shared_ptr<const Variables> varin_;
    };
  // -----------------------------------------------------------------------------

}  // namespace mom5cice5
#endif  // _MODEL_OBSFRACTION_H_
