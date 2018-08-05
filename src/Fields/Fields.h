/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_SRC_FIELDS_FIELDS_H_
#define SOCA_SRC_FIELDS_FIELDS_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "src/Geometry/Geometry.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "src/GetValuesTraj/GetValuesTraj.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
}

namespace ufo {
  class GeoVaLs;
}

namespace ioda {
  class Locations;
}

namespace soca {

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
      void getValues(const ioda::Locations &, const oops::Variables &,
                     ufo::GeoVaLs &) const;
      void getValues(const ioda::Locations &, const oops::Variables &,
                     ufo::GeoVaLs &, const GetValuesTraj &) const;
      void getValuesTL(const ioda::Locations &, const oops::Variables &,
                       ufo::GeoVaLs &, const GetValuesTraj &) const;
      void getValuesAD(const ioda::Locations &, const oops::Variables &,
                       const ufo::GeoVaLs &, const GetValuesTraj &);

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
      oops::Variables vars_;
      util::DateTime time_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_SRC_FIELDS_FIELDS_H_
