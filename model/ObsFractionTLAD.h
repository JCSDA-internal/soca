
#ifndef MOM5CICE5_MODEL_OBSFRACTIONTLAD_H_
#define MOM5CICE5_MODEL_OBSFRACTIONTLAD_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "model/LinearObsOp.h"
#include "model/ObsSpace.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace mom5cice5 {
  class Gom;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVec;

// -----------------------------------------------------------------------------
/// Sea-ice fraction observation for  model.
/*!
 *  ObsFractionTLAD for  model inherits from ObsEquivalent.
 */

class ObsFractionTLAD : public LinearObsOp, private util::ObjectCounter<ObsFractionTLAD> {
 public:
  static const std::string classname() {return "mom5cice5::ObsFractionTLAD";}

  ObsFractionTLAD(const ObsSpace &, const int &);
  virtual ~ObsFractionTLAD();

// Obs Operators
  void setTrajectory(const Gom &, const ObsBias &);
  void obsEquivTL(const Gom &, ObsVec &, const ObsBiasIncrement &) const;
  void obsEquivAD(Gom &, const ObsVec &, ObsBiasIncrement &) const;

// Other
  boost::shared_ptr<const Variables> variables() const {return varin_;}

  int& toFortran() {return keyOperStrm_;}
  const int& toFortran() const {return keyOperStrm_;}

 private:
  int keyOperStrm_;
  boost::shared_ptr<const Variables> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_OBSFRACTIONTLAD_H_
