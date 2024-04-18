/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_STATE_STATE_H_
#define SOCA_STATE_STATE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "soca/Fortran.h"
#include "soca/Fields/Fields.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace ufo {
  class GeoVaLs;
}
namespace soca {
  class Geometry;
  class Increment;
}

//-----------------------------------------------------------------------------

namespace soca {

  /// SOCA model state
  /*!
   * A State contains everything that is needed to propagate the state
   * forward in time.
   */
  class State : public Fields,
                private util::ObjectCounter<State> {
   public:
      static const std::string classname() {return "soca::State";}

      /// Constructor, destructor
      State(const Geometry &, const oops::Variables &,
            const util::DateTime &);
      State(const Geometry &, const eckit::Configuration &);
      State(const Geometry &, const State &);
      State(const oops::Variables &, const State &);
      State(const State &);
      virtual ~State();
      State & operator=(const State &);

      /// Rotations
      void rotate2north(const oops::Variables &, const oops::Variables &) const;
      void rotate2grid(const oops::Variables &, const oops::Variables &) const;

      /// Staggered grid interpolation
      void tohgrid(const oops::Variables &, const oops::Variables &) const;
      void tocgrid(const oops::Variables &, const oops::Variables &) const;

      /// Logarithmic and exponential transformations
      void logtrans(const oops::Variables &) const;
      void expontrans(const oops::Variables &) const;

      /// Interactions with Increment
      State & operator+=(const Increment &);

      /// I/O and diagnostics
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
      double norm() const;

      int & toFortran() {return keyFlds_;}
      const int & toFortran() const {return keyFlds_;}

      /// Update the fields in variable changes
      void updateFields(const oops::Variables &);

      /// Other
      void zero();
      void accumul(const double &, const State &);

      /// ATLAS Interfaces
      void toFieldSet(atlas::FieldSet &) const;
      void fromFieldSet(const atlas::FieldSet &);

   private:
      F90flds keyFlds_;
  };
// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_STATE_STATE_H_
