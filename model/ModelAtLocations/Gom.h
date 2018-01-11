
#ifndef SOCA_MODEL_GOM_H_
#define SOCA_MODEL_GOM_H_

#include <ostream>
#include <string>

#include "model/Fortran.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace oops {
  class Variables;
}

namespace soca {
  class ObsSpace;
  class Geometry;
  //class Variables;
  class Loc;
  
  /// Gom class to handle local model values for  model.

  class Gom : public util::Printable,
    private util::ObjectCounter<Gom> {
  public:
      static const std::string classname() {return "soca::Gom";}

      //Gom(const ObsSpace &, const oops::Variables &,
      //const util::DateTime &, const util::DateTime &);//, const Geometry &);
      Gom(const Loc &, const oops::Variables &);
      Gom(const eckit::Configuration &);
      
      explicit Gom(): keyGom_(0) {}
      explicit Gom(int & fgom): keyGom_(fgom) {}

      ~Gom();

      void zero();
      void random();
      Gom & operator*=(const double &);
      double dot_product_with(const Gom &) const;
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
  
      int & toFortran() {return keyGom_;}
      const int & toFortran() const {return keyGom_;}

  private:
      void print(std::ostream &) const;
      F90goms keyGom_;
    };

}  // namespace soca

#endif  // SOCA_MODEL_GOM_H_
