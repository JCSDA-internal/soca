
#ifndef SOCA_MODEL_SOCASTATE_H_
#define SOCA_MODEL_SOCASTATE_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "model/Fields/Fields.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

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
  //class Gom;
  //class Loc;
  class Geometry;
  class Increment;

  /// SOCA model state
  /*!
   * A State contains everything that is needed to propagate the state
   * forward in time.
   */

  // -----------------------------------------------------------------------------
  class State : public util::Printable,
    private util::ObjectCounter<State> {
  public:
      static const std::string classname() {return "soca::State";}

      /// Constructor, destructor
      State(const Geometry &, const oops::Variables &, const util::DateTime &);  // Is it used?
      State(const Geometry &, const eckit::Configuration &);
      State(const Geometry &, const State &);
      State(const State &);
      virtual ~State();
      State & operator=(const State &);

      /// Interpolate to observation location
      void interpolate(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;

      /// Interpolate full fields
      ///  void changeResolution(const State & xx);

      /// Interactions with Increment
      State & operator+=(const Increment &);

      /// I/O and diagnostics
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
      double norm() const {return fields_->norm();}
      const util::DateTime & validTime() const {return fields_->time();}
      util::DateTime & validTime() {return fields_->time();}

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

      void zero();
      void accumul(const double &, const State &);

  private:
      void print(std::ostream &) const;
      boost::scoped_ptr<Fields> fields_;
      boost::scoped_ptr<Fields> stash_;
    };
  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_SOCASTATE_H_
