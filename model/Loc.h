#ifndef MOM5CICE5_MODEL_LOC_H_
#define MOM5CICE5_MODEL_LOC_H_

#include <ostream>
#include <string>

#include "model/ObsSpace.h"
#include "model/Fortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace mom5cice5 {

  /// Loc class to handle locations for MOM5CICE5 model.

  class Loc : public util::Printable,
    private util::ObjectCounter<Loc> {
  public:
      static const std::string classname() {return "mom5cice5::Loc";}

      Loc(const ObsSpace & ot,
	  const util::DateTime & t1, const util::DateTime & t2) {
	keyLoc_ = ot.locations(t1, t2);
      }

      ~Loc() {mom5cice5_loc_delete_f90(keyLoc_);}

      int nobs() const {
	int nobs;
	mom5cice5_loc_nobs_f90(keyLoc_, nobs);
	return nobs;
      }

      int toFortran() const {return keyLoc_;}
  private:
      void print(std::ostream & os) const {
	os << "Loc::print not implemented";
      }
      int keyLoc_;
    };

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_LOC_H_
