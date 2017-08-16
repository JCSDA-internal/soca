
#ifndef MOM5CICE5_MODEL_MOM5CICE5GEOMETRY_H_
#define MOM5CICE5_MODEL_MOM5CICE5GEOMETRY_H_

#include <ostream>
#include <string>
#include <vector>

#include "model/Fortran.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace mom5cice5 {

  // -----------------------------------------------------------------------------
  /// Geometry handles geometry for MOM5CICE5 model.

  class Geometry : public util::Printable,
    private util::ObjectCounter<Geometry> {
  public:
      static const std::string classname() {return "mom5cice5::Geometry";}

      explicit Geometry(const eckit::Configuration &);
      Geometry(const Geometry &);
      ~Geometry();

      //std::vector<int> getDims() const;
      std::vector<double> getLats() const;
      std::vector<double> getLons() const;
      std::vector<double> getLevs() const;
      //std::vector<double> getArea() const;
      std::vector<int> getMask(const int &) const;

      int& toFortran() {return keyGeom_;}
      const int& toFortran() const {return keyGeom_;}

  private:
      Geometry & operator=(const Geometry &);
      void print(std::ostream &) const;
      int keyGeom_;
    };
  // -----------------------------------------------------------------------------

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_MOM5CICE5GEOMETRY_H_
