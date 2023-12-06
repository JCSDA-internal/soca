/*
 * (C) Copyright 2017-2023 UCAR
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

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & vt)
    : time_(vt), vars_(vars), geom_(geom)
  {
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_);
    Log::trace() << "State::State created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const eckit::Configuration & conf)
    : time_(),
      vars_(conf, "state variables"),
      geom_(geom)
  {
    util::DateTime * dtp = &time_;
    oops::Variables vars(vars_);
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars);

    if (conf.has("analytic init")) {
      std::string dt;
      conf.get("date", dt);
      time_ = util::DateTime(dt);
      soca_state_analytic_f90(toFortran(), &conf, &dtp);
    } else {
      soca_state_read_file_f90(toFortran(), &conf, &dtp);
    }
    Log::trace() << "State::State created and read in." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(geom)
  {
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_);
    soca_state_change_resol_f90(toFortran(), other.keyFlds_);
    Log::trace() << "State::State created by interpolation." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(other.geom_)
  {
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_);
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
  /// Staggered grid interpolation
  // -----------------------------------------------------------------------------
  void State::tohgrid(const oops::Variables & u,
                      const oops::Variables & v) const {
    Log::trace() << "State::State interpolate vector to h-grid."
                 << std::endl;
    soca_state_tohgrid_f90(toFortran());
  }
  // -----------------------------------------------------------------------------
  void State::tocgrid(const oops::Variables & u,
                      const oops::Variables & v) const {
    Log::trace() << "State::State interpolate vector to c-grid. NOT IMPLEMENTED"
                 << std::endl;
  }
  // -----------------------------------------------------------------------------
  /// Interactions with Increments
  // -----------------------------------------------------------------------------
  State & State::operator+=(const Increment & dx) {
    ASSERT(validTime() == dx.validTime());
    // Interpolate increment to analysis grid
    Increment dx_hr(geom_, dx);

    // Add increment to background state
    soca_state_add_incr_f90(toFortran(), dx_hr.toFortran());
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
    soca_state_sizes_f90(toFortran(), n0, n0, n0, nf);
    std::vector<double> zstat(3*nf);
    soca_state_gpnorm_f90(toFortran(), nf, zstat[0]);
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
  /// Serialization
  // -----------------------------------------------------------------------------
  size_t State::serialSize() const {
    // Field
    size_t nn;
    soca_state_serial_size_f90(toFortran(), geom_.toFortran(), nn);

    // Magic factor
    nn += 1;

    // Date and time
    nn += time_.serialSize();
    return nn;
  }
  // -----------------------------------------------------------------------------
  constexpr double SerializeCheckValue = -54321.98765;
  void State::serialize(std::vector<double> & vect) const {
    // Serialize the field
    size_t nn;
    soca_state_serial_size_f90(toFortran(), geom_.toFortran(), nn);
    std::vector<double> vect_field(nn, 0);
    vect.reserve(vect.size() + nn + 1 + time_.serialSize());
    soca_state_serialize_f90(toFortran(), geom_.toFortran(), nn,
                             vect_field.data());
    vect.insert(vect.end(), vect_field.begin(), vect_field.end());

    // Magic value placed in serialization; used to validate deserialization
    vect.push_back(SerializeCheckValue);

    // Serialize the date and time
    time_.serialize(vect);
  }
  // -----------------------------------------------------------------------------
  void State::deserialize(const std::vector<double> & vect, size_t & index) {
    // Deserialize the field
    soca_state_deserialize_f90(toFortran(), geom_.toFortran(), vect.size(),
                               vect.data(), index);

    // Use magic value to validate deserialization
    ASSERT(vect.at(index) == SerializeCheckValue);
    ++index;

    // Deserialize the date and time
    time_.deserialize(vect, index);
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
  void State::updateFields(const oops::Variables & vars) {
    // Update local variables
    vars_ = vars;
    // Update field data
    soca_state_update_fields_f90(toFortran(), vars_);
  }
  // -----------------------------------------------------------------------------
  /// Logarithmic and exponential transformations
  // -----------------------------------------------------------------------------
  void State::logtrans(const oops::Variables & trvar) const {
    Log::trace() << "State::State apply logarithmic transformation."
                 << std::endl;
    soca_state_logtrans_f90(toFortran(), trvar);
  }
  // -----------------------------------------------------------------------------
  void State::expontrans(const oops::Variables & trvar) const {
    Log::trace() << "State::State apply exponential transformation."
    << std::endl;
    soca_state_expontrans_f90(toFortran(), trvar);
  }

  // -----------------------------------------------------------------------------
  const util::DateTime & State::validTime() const {return time_;}
  // -----------------------------------------------------------------------------
  util::DateTime & State::validTime() {return time_;}
  // -----------------------------------------------------------------------------
  const Geometry & State::geometry() const {return geom_;}
  // -----------------------------------------------------------------------------

  void State::toFieldSet(atlas::FieldSet &fset) const {
    soca_state_to_fieldset_f90(toFortran(), vars_, fset.get());
    fset.haloExchange();
  }

  // -----------------------------------------------------------------------------

  void State::fromFieldSet(const atlas::FieldSet &fs) {
    soca_state_from_fieldset_f90(toFortran(), vars_, fs.get());
  }
}  // namespace soca
