/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iomanip>
#include <vector>

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"
#include "soca/State/StateFortran.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & vt)
    : time_(vt), vars_(vars), geom_(new Geometry(geom))
  {
    soca_state_create_f90(keyFlds_, geom_->toFortran(), vars_);
    Log::trace() << "State::State created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars,
               const eckit::Configuration & file)
    : time_(), vars_(vars), geom_(new Geometry(geom))
  {
    util::DateTime * dtp = &time_;
    soca_state_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_state_read_file_f90(toFortran(), &file, &dtp);
    Log::trace() << "State::State created and read in." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(new Geometry(geom))
  {
    soca_state_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_state_change_resol_f90(toFortran(), other.keyFlds_);
    Log::trace() << "State::State created by interpolation." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(new Geometry(*other.geom_))
  {
    soca_state_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_state_copy_f90(toFortran(), other.toFortran());
    Log::trace() << "State::State copied." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::~State() {
    soca_state_delete_f90(toFortran());
    Log::trace() << "State::State destructed." << std::endl;
  }
  // -----------------------------------------------------------------------------
  /// Basic operators
  // -----------------------------------------------------------------------------
  State & State::operator=(const State & rhs) {
    time_ = rhs.time_;
    soca_state_copy_f90(toFortran(), rhs.toFortran());
    return *this;
  }

  // -----------------------------------------------------------------------------
  /// Rotations
  // -----------------------------------------------------------------------------
  void State::rotate2north(const oops::Variables & u,
                           const oops::Variables & v) const {
    Log::trace() << "State::State rotate from logical to geographical North."
                 << std::endl;
    soca_state_rotate2north_f90(toFortran(), u, v);
  }
  // -----------------------------------------------------------------------------
  void State::rotate2grid(const oops::Variables & u,
                          const oops::Variables & v) const {
    Log::trace() << "State::State rotate from geographical to logical North."
    << std::endl;
    soca_state_rotate2grid_f90(toFortran(), u, v);
  }
  // -----------------------------------------------------------------------------
  /// Interactions with Increments
  // -----------------------------------------------------------------------------
  State & State::operator+=(const Increment & dx) {
    ASSERT(validTime() == dx.validTime());
    soca_state_add_incr_f90(toFortran(), dx.toFortran());
    return *this;
  }
  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------
  void State::read(const eckit::Configuration & files) {
    Log::trace() << "State::State read started." << std::endl;
    util::DateTime * dtp = &time_;
    soca_state_read_file_f90(toFortran(), &files, &dtp);
    Log::trace() << "State::State read done." << std::endl;
  }
  // -----------------------------------------------------------------------------
  void State::write(const eckit::Configuration & files) const {
    const util::DateTime * dtp = &time_;
    soca_state_write_file_f90(toFortran(), &files, &dtp);
  }
  // -----------------------------------------------------------------------------
  void State::print(std::ostream & os) const {
    os << std::endl << "  Valid time: " << validTime();
    int n0, nf;
    soca_state_sizes_f90(toFortran(), n0, n0, n0, n0, n0, nf);
    std::vector<double> zstat(3*nf);
    soca_state_gpnorm_f90(toFortran(), nf, zstat[0]);
    for (int jj = 0; jj < nf; ++jj) {
      // TODO(travis) remove this once answers ready to be changed
      if (vars_[jj] == "mld" || vars_[jj] == "layer_depth") {
        continue;
      }

      os << std::endl << std::right << std::setw(7) << vars_[jj]
         << "   min="  <<  std::fixed << std::setw(12) <<
                           std::right << zstat[3*jj]
         << "   max="  <<  std::fixed << std::setw(12) <<
                           std::right << zstat[3*jj+1]
         << "   mean=" <<  std::fixed << std::setw(12) <<
                           std::right << zstat[3*jj+2];
    }
  }

  // -----------------------------------------------------------------------------
  /// For accumulator
  // -----------------------------------------------------------------------------
  void State::zero() {
    soca_state_zero_f90(toFortran());
  }
  // -----------------------------------------------------------------------------
  void State::accumul(const double & zz, const State & xx) {
    soca_state_axpy_f90(toFortran(), zz, xx.toFortran());
  }
  // -----------------------------------------------------------------------------
  double State::norm() const {
    double zz = 0.0;
    soca_state_rms_f90(toFortran(), zz);
    return zz;
  }
  // -----------------------------------------------------------------------------
  const util::DateTime & State::validTime() const {return time_;}
  // -----------------------------------------------------------------------------
  util::DateTime & State::validTime() {return time_;}
  // -----------------------------------------------------------------------------
  boost::shared_ptr<const Geometry> State::geometry() const {return geom_;}
  // -----------------------------------------------------------------------------
}  // namespace soca
