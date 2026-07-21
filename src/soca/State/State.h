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

#include "eckit/mpi/Comm.h"

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
      void transpose(const State & /*DistState*/, const eckit::mpi::Comm & /*global*/,
         const int /*ensNum*/, const int /*transNum*/) {
         throw eckit::NotImplemented("Soca State::transpose not implemented", Here());
      }

      /// Rotations
      void rotate2north(const oops::Variables &, const oops::Variables &);
      void rotate2grid(const oops::Variables &, const oops::Variables &);

      /// Staggered grid interpolation
      void tohgrid(const oops::Variables &, const oops::Variables &);
      void tocgrid(const oops::Variables &, const oops::Variables &);

      /// Logarithmic and exponential transformations
      void logtrans(const oops::Variables &);
      void expontrans(const oops::Variables &);

      /// Interactions with Increment
      State & operator+=(const Increment &);

      /// I/O and diagnostics
      void read(const eckit::Configuration &);
      void write(const eckit::Configuration &) const;

      /// Bulk parallel write across ensemble members. Each member is gathered
      /// to a strided writer PE (assigned by soca_io_ensemble_root_pe in
      /// soca_io_mod) and written in phase 2 with no MPI, so per-member
      /// netCDF writes happen concurrently across writer PEs. See
      /// soca_io_writers_commit_ensemble for the staging.
      static void writeEnsemble(const std::vector<const State*> &,
                                const std::vector<eckit::LocalConfiguration> &);

      /// Bulk read across ensemble members via soca_io_readers_commit_ensemble.
      /// Members are read in parallel (rotated reader PEs), honoring the
      /// geometry.io.'ensemble read' (scatter vs strided) and 'async mpi'
      /// toggles. See soca_io_readers_commit_ensemble for the staging.
      static void readEnsemble(const std::vector<State*> &,
                               const std::vector<eckit::LocalConfiguration> &);

      int & toFortran() {return keyFlds_;}
      const int & toFortran() const {return keyFlds_;}

      /// Update the fields in variable changes
      void updateFields(const oops::Variables &);

   private:
      F90flds keyFlds_;
  };
// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_STATE_STATE_H_
