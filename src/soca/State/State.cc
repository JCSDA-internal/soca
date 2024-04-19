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
#include "soca/VariableChange/VariableChange.h"

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"

#include "ufo/GeoVaLs.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & vt)
    : Fields(geom, vars, vt)
  {
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    syncToFieldset();
    Log::trace() << "State::State created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const eckit::Configuration & conf)
    : Fields(geom, oops::Variables(conf, "state variables"), util::DateTime())
  {
    util::DateTime * dtp = &time_;
    oops::Variables vars(vars_);
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars, fieldSet_.get());
    syncToFieldset();

    if (conf.has("analytic init")) {
      std::string dt;
      conf.get("date", dt);
      time_ = util::DateTime(dt);
      soca_state_analytic_f90(toFortran(), &conf, &dtp);
    } else {
      soca_state_read_file_f90(toFortran(), &conf, &dtp);
    }
    syncToFieldset();
    Log::trace() << "State::State created and read in." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const State & other)
    : Fields(geom, other.vars_, other.time_)
  {
    other.syncFromFieldset();
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    soca_state_change_resol_f90(toFortran(), other.keyFlds_);
    syncToFieldset();
    Log::trace() << "State::State created by interpolation." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const oops::Variables & vars, const State & other) : State(other)
  {
    // TODO(Travis, maybe) The variable change needs to go here
    //  (U/V rotate, etc), but since we don't really use it right now, this
    //  can wait until after the variable changes get cleaned up, (and
    //  after we finally implement the model naming convention??)

    // eckit::LocalConfiguration varChangeConfig;
    // varChangeConfig.set("variable change name", "Model2Ana");
    // VariableChange model2ana(varChangeConfig, geom_);
    // model2ana.changeVar(*this, vars);
    Log::trace() << "State::State created with variable change." << std::endl;
  }
  // -----------------------------------------------------------------------------
  State::State(const State & other)
    : Fields(other.geom_, other.vars_, other.time_)
  {
    other.syncFromFieldset();
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    soca_state_copy_f90(toFortran(), other.toFortran());
    syncToFieldset();
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
    rhs.syncFromFieldset();
    syncFromFieldset();
    soca_state_copy_f90(toFortran(), rhs.toFortran());
    syncToFieldset();
    return *this;
  }
  // -----------------------------------------------------------------------------
  /// Rotations
  // -----------------------------------------------------------------------------
  void State::rotate2north(const oops::Variables & u,
                           const oops::Variables & v) {
    Log::trace() << "State::State rotate from logical to geographical North."
                 << std::endl;
    soca_state_rotate2north_f90(toFortran(), u, v);
    syncToFieldset();
  }
  // -----------------------------------------------------------------------------
  void State::rotate2grid(const oops::Variables & u,
                          const oops::Variables & v) {
    Log::trace() << "State::State rotate from geographical to logical North."
    << std::endl;
    soca_state_rotate2grid_f90(toFortran(), u, v);
    syncToFieldset();
  }
  // -----------------------------------------------------------------------------
  /// Staggered grid interpolation
  // -----------------------------------------------------------------------------
  void State::tohgrid(const oops::Variables & u,
                      const oops::Variables & v) {
    Log::trace() << "State::State interpolate vector to h-grid."
                 << std::endl;
    soca_state_tohgrid_f90(toFortran());
    syncToFieldset();
  }
  // -----------------------------------------------------------------------------
  void State::tocgrid(const oops::Variables & u,
                      const oops::Variables & v) {
    Log::trace() << "State::State interpolate vector to c-grid. NOT IMPLEMENTED"
                 << std::endl;
    syncToFieldset();
  }
  // -----------------------------------------------------------------------------
  /// Interactions with Increments
  // -----------------------------------------------------------------------------
  State & State::operator+=(const Increment & dx) {
    ASSERT(validTime() == dx.validTime());
    // Interpolate increment to analysis grid
    Increment dx_hr(geom_, dx);

    // Add increment to background state
    atlas::FieldSet fs2; dx.toFieldSet(fs2);
    util::addFieldSets(fieldSet_, fs2);
    syncFromFieldset();
    return *this;
  }
  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------
  void State::read(const eckit::Configuration & files) {
    Log::trace() << "State::State read started." << std::endl;
    util::DateTime * dtp = &time_;
    soca_state_read_file_f90(toFortran(), &files, &dtp);
    syncToFieldset();
    Log::trace() << "State::State read done." << std::endl;
  }
  // -----------------------------------------------------------------------------
  void State::write(const eckit::Configuration & files) const {
    const util::DateTime * dtp = &time_;
    syncFromFieldset();
    soca_state_write_file_f90(toFortran(), &files, &dtp);
  }

  // -----------------------------------------------------------------------------
  /// For accumulator
  // -----------------------------------------------------------------------------
  void State::zero() {
    util::zeroFieldSet(fieldSet_);
    syncFromFieldset();
  }
  // -----------------------------------------------------------------------------
  void State::accumul(const double & zz, const State & xx) {
    atlas::FieldSet fs1, fs2; xx.toFieldSet(fs1);
    fs2 = util::copyFieldSet(fs1);
    util::multiplyFieldSet(fs2, zz);
    util::addFieldSets(fieldSet_, fs2);
    syncFromFieldset();
  }
  // -----------------------------------------------------------------------------
  void State::updateFields(const oops::Variables & vars) {
    // Update local variables
    vars_ = vars;
    // Update field data
    syncFromFieldset();
    soca_state_update_fields_f90(toFortran(), vars_);
    syncToFieldset();
  }
  // -----------------------------------------------------------------------------
  /// Logarithmic and exponential transformations
  // -----------------------------------------------------------------------------
  void State::logtrans(const oops::Variables & trvar) {
    Log::trace() << "State::State apply logarithmic transformation."
                 << std::endl;
    syncFromFieldset();
    soca_state_logtrans_f90(toFortran(), trvar);
    syncToFieldset();
  }
  // -----------------------------------------------------------------------------
  void State::expontrans(const oops::Variables & trvar) {
    Log::trace() << "State::State apply exponential transformation."
    << std::endl;
    syncFromFieldset();
    soca_state_expontrans_f90(toFortran(), trvar);
    syncToFieldset();
  }

  // -----------------------------------------------------------------------------

  void State::toFieldSet(atlas::FieldSet &fset) const {
    util::copyFieldSet(fieldSet_, fset);
  }

  // -----------------------------------------------------------------------------

  void State::fromFieldSet(const atlas::FieldSet &fs) {
    util::copyFieldSet(fs, fieldSet_);
    soca_state_from_fieldset_f90(toFortran(), vars_, fs.get());
  }

  void State::syncFromFieldset() const{
    soca_state_from_fieldset_f90(toFortran(), vars_, fieldSet_.get());
  }
  void State::syncToFieldset() const{
    soca_state_to_fieldset_f90(toFortran(), vars_, fieldSet_.get());
  }

}  // namespace soca
