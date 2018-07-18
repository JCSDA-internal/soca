/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SOCA_MODEL_SOCAINCREMENT_H_
#define SOCA_MODEL_SOCAINCREMENT_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "src/Fields/Fields.h"
#include "src/GetValuesTraj/GetValuesTraj.h"
#include "src/Geometry/Geometry.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/dot_product.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
  class Variables;
}

namespace ufo {
  class GeoVaLs;
}

namespace ioda {
  class Locations;
}

namespace soca {
  class ModelBiasIncrement;
  class ErrorCovariance;  
  class State;

  /// Increment Class: Difference between two states
  /*!
   *  Some fields that are present in a State may not be present in
   *  an Increment. The Increment contains everything that is needed by
   *  the tangent-linear and adjoint models.
   */
  
  // -----------------------------------------------------------------------------

  class Increment : public oops::GeneralizedDepartures,
    public util::Printable,
    private util::ObjectCounter<Increment> {
  public:
      static const std::string classname() {return "soca::Increment";}

      /// Constructor, destructor
      Increment(const Geometry &, const oops::Variables &, const util::DateTime &);
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
      
      /// Interpolate to observation location
      void getValuesTL(const ioda::Locations &, const oops::Variables &,ufo::GeoVaLs &, const GetValuesTraj &) const;
      void getValuesAD(const ioda::Locations &, const oops::Variables &, const ufo::GeoVaLs &, const GetValuesTraj &);

      /// I/O and diagnostics
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
      double norm() const {return fields_->norm();}
      const util::DateTime & validTime() const {return fields_->time();}
      util::DateTime & validTime() {return fields_->time();}
      void updateTime(const util::Duration & dt) {fields_->time() += dt;}


      /// Define and convert to/from unstructured grid
      void define(oops::UnstructuredGrid &) const;      
      void convert_to(oops::UnstructuredGrid &) const;
      void convert_from(const oops::UnstructuredGrid &);
  
      /// Access to fields
      Fields & fields() {return *fields_;}
      const Fields & fields() const {return *fields_;}

      boost::shared_ptr<const Geometry> geometry() const {
	return fields_->geometry();
      }

      /// Other
      void activateModel();
      void deactivateModel();

      void accumul(const double &, const State &);

      /// Data
  private:
      void print(std::ostream &) const;
      boost::scoped_ptr<Fields> fields_;
      boost::scoped_ptr<Fields> stash_;
    };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_SOCAINCREMENT_H_
