/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>

#include "soca/Fields/FieldsFortran.h"
#include "soca/Geometry/Geometry.h"
#include "soca/GetValuesTraj/GetValuesTraj.h"
#include "soca/Increment/Increment.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
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
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    Log::trace() << "State::State created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars,
               const eckit::Configuration & file)
    : time_(), vars_(vars), geom_(new Geometry(geom))
  {
    util::DateTime * dtp = &time_;
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_field_read_file_f90(toFortran(), &file, &dtp);
    Log::trace() << "State::State created and read in." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(new Geometry(geom))
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_field_copy_f90(toFortran(), other.keyFlds_);
    Log::trace() << "State::State created by interpolation." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(new Geometry(*other.geom_))
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_field_copy_f90(toFortran(), other.toFortran());
    Log::trace() << "State::State copied." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::~State() {
    soca_field_delete_f90(toFortran());
    Log::trace() << "State::State destructed." << std::endl;
  }
  // -----------------------------------------------------------------------------
  /// Basic operators
  // -----------------------------------------------------------------------------
  State & State::operator=(const State & rhs) {
    time_ = rhs.time_;
    soca_field_copy_f90(toFortran(), rhs.toFortran());
    return *this;
  }
  // -----------------------------------------------------------------------------
  /// Interpolate to observation location
  // -----------------------------------------------------------------------------
  void State::getValues(const ufo::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & cols) const {
    if (geom_->getAtmInit())
      {
        // Get atm geovals
        // The variables in vars that are also defined in soca will be
        // over-written in the interpolation call bellow
        getValuesFromFile(locs, vars, cols);
        }
    // Get ocean geovals
    soca_field_interp_nl_f90(toFortran(), locs.toFortran(), vars,
                             cols.toFortran());
  }
  // -----------------------------------------------------------------------------
  void State::getValues(const ufo::Locations & locs,
                        const oops::Variables & vars,
                        ufo::GeoVaLs & cols,
                        GetValuesTraj & traj) const {
    soca_field_interp_nl_traj_f90(toFortran(), locs.toFortran(), vars,
                                  cols.toFortran(), traj.toFortran());
  }
  // -----------------------------------------------------------------------------
  /// Read Interpolated GeoVaLs from file
  // -----------------------------------------------------------------------------
  void State::getValuesFromFile(const ufo::Locations & locs,
                                const oops::Variables & vars,
                                ufo::GeoVaLs & atmgom) const {
    // Get Atm configuration
    eckit::LocalConfiguration conf(geom_->getAtmConf());

    // Get Time Bounds
    util::DateTime bgn = util::DateTime(conf.getString("notocean.date_begin"));
    util::DateTime end = util::DateTime(conf.getString("notocean.date_end"));

    // Create the Atmospheric Geometry in Observation Space
    eckit::LocalConfiguration confatmobs(conf, "notocean.ObsSpace");
    ioda::ObsSpace atmobs(confatmobs, geom_->getComm(),
                                bgn, end);

    // Get GeoVaLs from file
    eckit::LocalConfiguration confatm(conf, "notocean");
    atmgom.read(confatm, atmobs);
  }

  // -----------------------------------------------------------------------------
  /// Rotations
  // -----------------------------------------------------------------------------
  void State::rotate2north(const oops::Variables & u,
                           const oops::Variables & v) const {
    Log::trace() << "State::State rotate from logical to geographical North."
                 << std::endl;
    fields_->rotate2north(u, v);
  }
  // -----------------------------------------------------------------------------
  void State::rotate2grid(const oops::Variables & u,
                          const oops::Variables & v) const {
    Log::trace() << "State::State rotate from geographical to logical North."
    << std::endl;
    fields_->rotate2grid(u, v);
  }
  // -----------------------------------------------------------------------------
  /// Interactions with Increments
  // -----------------------------------------------------------------------------
  State & State::operator+=(const Increment & dx) {
    ASSERT(validTime() == dx.validTime());
    soca_field_add_incr_f90(toFortran(), dx.toFortran());
    return *this;
  }
  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------
  void State::read(const eckit::Configuration & files) {
    Log::trace() << "State::State read started." << std::endl;
    util::DateTime * dtp = &time_;
    soca_field_read_file_f90(toFortran(), &files, &dtp);
    Log::trace() << "State::State read done." << std::endl;
  }
  // -----------------------------------------------------------------------------
  void State::write(const eckit::Configuration & files) const {
    const util::DateTime * dtp = &time_;
    soca_field_write_file_f90(toFortran(), &files, &dtp);
  }
  // -----------------------------------------------------------------------------
  void State::print(std::ostream & os) const {
    os << std::endl << "  Valid time: " << validTime();
    int n0, nf;
    soca_field_sizes_f90(toFortran(), n0, n0, n0, n0, n0, nf);
    std::vector<double> zstat(3*nf);
    soca_field_gpnorm_f90(toFortran(), nf, zstat[0]);
    for (int jj = 0; jj < nf; ++jj) {
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
    soca_field_zero_f90(toFortran());
  }
  // -----------------------------------------------------------------------------
  void State::accumul(const double & zz, const State & xx) {
    soca_field_axpy_f90(toFortran(), zz, xx.toFortran());
  }
  // -----------------------------------------------------------------------------
  double State::norm() const {
    double zz = 0.0;
    soca_field_rms_f90(toFortran(), zz);
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
