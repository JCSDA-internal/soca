
#ifndef MOM5CICE5_MODEL_OBSFRACTIONTLAD_H_
#define MOM5CICE5_MODEL_OBSFRACTIONTLAD_H_

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "model/Traits.h"
#include "oops/interface/LinearObsOperBase.h"
#include "model/ObsSpace/ObsSpace.h"
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

  class ObsFractionTLAD : public oops::LinearObsOperBase<Traits>, private util::ObjectCounter<ObsFractionTLAD> {
  public:
    static const std::string classname() {return "mom5cice5::ObsFractionTLAD";}

    //ObsFractionTLAD(const ObsSpace &, const int &);
    ObsFractionTLAD(const ObsSpace &, const eckit::Configuration &);    
    virtual ~ObsFractionTLAD();

    // Obs Operators
    void setTrajectory(const Gom &, const ObsBias &);
    void obsEquivTL(const Gom &, ObsVec &, const ObsBiasIncrement &) const override;
    void obsEquivAD(Gom &, const ObsVec &, ObsBiasIncrement &) const override;

    // Other
    boost::shared_ptr<const Variables> variables() const override {return varin_;}

    int& toFortran() {return keyOperFraction_;}
    const int& toFortran() const {return keyOperFraction_;}

  private:
    void print(std::ostream &) const override;    
    F90hop keyOperFraction_;
    boost::shared_ptr<const Variables> varin_;
  };
  // -----------------------------------------------------------------------------

}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_OBSFRACTIONTLAD_H_
