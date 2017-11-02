
#ifndef MOM5CICE5_MODEL_GOM_H_
#define MOM5CICE5_MODEL_GOM_H_

#include <ostream>
#include <string>

#include "model/Fortran.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace mom5cice5 {
  class ObsSpace;
  class Geometry;
  class Variables;

  /// Gom class to handle local model values for  model.

  class Gom : public util::Printable,
    private util::ObjectCounter<Gom> {
  public:
      static const std::string classname() {return "mom5cice5::Gom";}

      Gom(const ObsSpace &, const Variables &,
	  const util::DateTime &, const util::DateTime &, const Geometry &);

      explicit Gom(): keyGom_(0) {}
      explicit Gom(int & fgom): keyGom_(fgom) {}

      ~Gom();

      void zero();
      double dot_product_with(const Gom & other) const;

      int & toFortran() {return keyGom_;}
      const int & toFortran() const {return keyGom_;}

  private:
      void print(std::ostream &) const;
      int keyGom_;
    };

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_GOM_H_
