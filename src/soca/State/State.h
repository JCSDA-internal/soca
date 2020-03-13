/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_STATE_STATE_H_
#define SOCA_STATE_STATE_H_

#include <memory>
#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace oops {
  class Variables;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace soca {
  class Fields;
  class Geometry;
  class GetValuesTraj;
  class Increment;
}

//-----------------------------------------------------------------------------

namespace soca {

  /// SOCA model state
  /*!
   * A State contains everything that is needed to propagate the state
   * forward in time.
   */
  class State : public util::Printable,
    private util::ObjectCounter<State> {
   public:
      static const std::string classname() {return "soca::State";}

      /// Constructor, destructor
      State(const Geometry &, const oops::Variables &,
            const util::DateTime &);
      State(const Geometry &, const oops::Variables &,
            const eckit::Configuration &);
      State(const Geometry &, const State &);
      State(const State &);
      virtual ~State();
      State & operator=(const State &);

      /// Interpolate to observation location
      void getValues(const ufo::Locations &,
                     const oops::Variables &,
                     ufo::GeoVaLs &) const;
      void getValues(const ufo::Locations &,
                     const oops::Variables &,
                     ufo::GeoVaLs &,
                     GetValuesTraj &) const;

      /// Read interpolated GeoVaLs at observation location
      void getValuesFromFile(const ufo::Locations &,
                             const oops::Variables &,
                             ufo::GeoVaLs &) const;

      /// Grid 2 grid and rotations
      void rotate2north() const;
      void rotate2grid() const;
      void uvgrid2agrid() const;
      void agrid2uvgrid() const;

      /// Interpolate full fields
      ///  void changeResolution(const State & xx);

      /// Interactions with Increment
      State & operator+=(const Increment &);

      /// I/O and diagnostics
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
      double norm() const;
      const util::DateTime & validTime() const;
      util::DateTime & validTime();

      /// Access to fields
      Fields & fields();
      const Fields & fields() const;

      boost::shared_ptr<const Geometry> geometry() const;

      /// Other

      void zero();
      void accumul(const double &, const State &);

   private:
      void print(std::ostream &) const;
      std::unique_ptr<Fields> fields_;
      std::unique_ptr<Fields> stash_;
      boost::shared_ptr<const Geometry> geom_;
      oops::Variables vars_;
      util::DateTime time_;
  };
// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_STATE_STATE_H_
