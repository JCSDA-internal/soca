/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SOCA_INCREMENT_INCREMENT_H_
#define SOCA_INCREMENT_INCREMENT_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "soca/Fortran.h"

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class GridPoint;
  class UnstructuredGrid;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace soca {
  class Geometry;
  class GeometryIterator;
  class GetValuesTraj;
  class State;
}

// -----------------------------------------------------------------------------

namespace soca {

  /// Increment Class: Difference between two states
  /*!
   *  Some fields that are present in a State may not be present in
   *  an Increment. The Increment contains everything that is needed by
   *  the tangent-linear and adjoint models.
   */

  class Increment : public oops::GeneralizedDepartures,
    public util::Printable,
    private util::ObjectCounter<Increment> {
   public:
      static const std::string classname() {return "soca::Increment";}

      /// Constructor, destructor
      Increment(const Geometry &, const oops::Variables &,
                const util::DateTime &);
      Increment(const Geometry &, const Increment &);
      Increment(const Increment &, const bool);
      Increment(const Increment &);
      virtual ~Increment();

      /// Basic operators
      void diff(const State &, const State &);
      void zero();
      void zero(const util::DateTime &);
      Increment & operator =(const Increment &);
      Increment & operator+=(const Increment &);
      Increment & operator-=(const Increment &);
      Increment & operator*=(const double &);
      void axpy(const double &, const Increment &, const bool check = true);
      double dot_product_with(const Increment &) const;
      void schur_product_with(const Increment &);
      void random();
      void dirac(const eckit::Configuration &);

      /// Getpoint/Setpoint
      oops::GridPoint getPoint(const GeometryIterator &) const;
      void setPoint(const oops::GridPoint &, const GeometryIterator &);

      /// Interpolate to observation location
      void getValuesTL(const ufo::Locations &,
                       const oops::Variables &,
                       ufo::GeoVaLs &,
                       const GetValuesTraj &) const;
      void getValuesAD(const ufo::Locations &,
                       const oops::Variables &,
                       const ufo::GeoVaLs &,
                       const GetValuesTraj &);

      /// I/O and diagnostics
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
      double norm() const;
      const util::DateTime & validTime() const;
      util::DateTime & validTime();
      void updateTime(const util::Duration & dt);


      /// Unstructured grid
      void ug_coord(oops::UnstructuredGrid &) const;
      void field_to_ug(oops::UnstructuredGrid &, const int &) const;
      void field_from_ug(const oops::UnstructuredGrid &, const int &);

      /// Other
      void accumul(const double &, const State &);
      int & toFortran() {return keyFlds_;}
      const int & toFortran() const {return keyFlds_;}
      boost::shared_ptr<const Geometry> geometry() const;


      /// Data
   private:
      void print(std::ostream &) const;

      F90flds keyFlds_;
      oops::Variables vars_;
      util::DateTime time_;
      boost::shared_ptr<const Geometry> geom_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_INCREMENT_INCREMENT_H_
