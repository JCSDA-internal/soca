
#ifndef SOCA_MODEL_SOCAGEOMETRY_H_
#define SOCA_MODEL_SOCAGEOMETRY_H_

#include <ostream>
#include <string>
#include <vector>

#include "model/Fortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace soca {

  // -----------------------------------------------------------------------------
  /// Geometry handles geometry for SOCA model.

  class Geometry : public util::Printable,
    private util::ObjectCounter<Geometry> {
  public:
      static const std::string classname() {return "soca::Geometry";}

      explicit Geometry(const eckit::Configuration &);
      Geometry(const Geometry &);
      ~Geometry();

      int& toFortran() {return keyGeom_;}
      const int& toFortran() const {return keyGeom_;}

  private:
      Geometry & operator=(const Geometry &);
      void print(std::ostream &) const;
      int keyGeom_;
    };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_SOCAGEOMETRY_H_
