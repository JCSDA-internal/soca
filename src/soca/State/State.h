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
  class Locations;
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
  class State : public util::Printable,
                public util::Serializable,
    private util::ObjectCounter<State> {
   public:
      static const std::string classname() {return "soca::State";}

      /// Constructor, destructor
      State(const Geometry &, const oops::Variables &,
            const util::DateTime &);
      State(const Geometry &, const eckit::Configuration &);
      State(const Geometry &, const State &);
      State(const State &);
      virtual ~State();
      State & operator=(const State &);

      /// Needed by PseudoModel
      void updateTime(const util::Duration & dt) {time_ += dt;}

      /// Rotations
      void rotate2north(const oops::Variables &, const oops::Variables &) const;
      void rotate2grid(const oops::Variables &, const oops::Variables &) const;

      /// Logarithmic and exponential transformations
      void logtrans(const oops::Variables &) const;
      void expontrans(const oops::Variables &) const;

      /// Interactions with Increment
      State & operator+=(const Increment &);

      /// I/O and diagnostics
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;
      double norm() const;
      const util::DateTime & validTime() const;
      util::DateTime & validTime();

      /// Serialize and deserialize
      size_t serialSize() const override;
      void serialize(std::vector<double> &) const override;
      void deserialize(const std::vector<double> &, size_t &) override;


      int & toFortran() {return keyFlds_;}
      const int & toFortran() const {return keyFlds_;}
      std::shared_ptr<const Geometry> geometry() const;
      const oops::Variables & variables() const {return vars_;}
      const util::DateTime & time() const {return time_;}

      /// Update the fields in variable changes
      const bool hasFields(const oops::Variables &);
      void updateFields(const oops::Variables &);

      /// Other
      void zero();
      void accumul(const double &, const State &);

   private:
      void print(std::ostream &) const override;

      F90flds keyFlds_;

      std::shared_ptr<const Geometry> geom_;
      oops::Variables vars_;
      util::DateTime time_;
  };
// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_STATE_STATE_H_
