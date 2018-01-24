
#ifndef SOCA_MODEL_SOCAFIELDS_H_
#define SOCA_MODEL_SOCAFIELDS_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "model/Geometry/Geometry.h"
#include "model/Variables/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
  class Variables;  
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace soca {
  //class Loc;
  //class Gom;
  
  // -----------------------------------------------------------------------------
  /// Class to represent a FieldSet for the SOCA model
  class Fields : public util::Printable,
    private util::ObjectCounter<Fields> {
  public:
      static const std::string classname() {return "soca::Fields";}

      // Constructors and basic operators
      Fields(const Geometry &, const oops::Variables &, const util::DateTime &);
      Fields(const Fields &, const Geometry &);
      Fields(const Fields &, const oops::Variables &);
      Fields(const Fields &, const bool);
      Fields(const Fields &);
      ~Fields();

      void zero();
      void zero(const util::DateTime &);
      void dirac(const eckit::Configuration &);  
      Fields & operator=(const Fields &);
      Fields & operator+=(const Fields &);
      Fields & operator-=(const Fields &);
      Fields & operator*=(const double &);
      void axpy(const double &, const Fields &);
      double dot_product_with(const Fields &) const;
      void schur_product_with(const Fields &);
      void random();

      // Interpolate to given location
      void interpolate(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;
      void interpolateTL(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;
      void interpolateAD(const ufo::Locations &, const oops::Variables &, const ufo::GeoVaLs &);      
      //void interpolateTL(const Loc &, Gom &) const;
      //void interpolateAD(const Loc &, const Gom &);

      // Interpolate full fields
      void changeResolution(const Fields &);
      void add(const Fields &);
      void diff(const Fields &, const Fields &);
  
      // Convert to/from unstructured grid
      void convert_to(oops::UnstructuredGrid &) const;
      void convert_from(const oops::UnstructuredGrid &);
  
      // Utilities
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
      double norm() const;
      boost::shared_ptr<const Geometry> geometry() const {return geom_;}

      const util::DateTime & time() const {return time_;}
      util::DateTime & time() {return time_;}

      int & toFortran() {return keyFlds_;}
      const int & toFortran() const {return keyFlds_;}

      bool isForModel(const bool) const;

  private:
      void print(std::ostream &) const;
      F90flds keyFlds_;
      boost::shared_ptr<const Geometry> geom_;
      const Variables vars_;      
      util::DateTime time_;
    };
  // -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_MODEL_SOCAFIELDS_H_
