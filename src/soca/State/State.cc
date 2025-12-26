/*
 * (C) Copyright 2017-2024 UCAR
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
#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/GlobalInterpolator.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars, const util::DateTime & vt)
    : Fields(geom, vars, vt)
  {
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    Log::trace() << "State::State created." << std::endl;
  }

  // -----------------------------------------------------------------------------

  State::State(const Geometry & geom, const eckit::Configuration & conf)
    : State(geom, oops::Variables(conf, "state variables"), util::DateTime())
  {
    util::DateTime * dtp = &time_;

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
  // Resolution change
  State::State(const Geometry & geom, const State & other)
    : Fields(geom, other.vars_, other.time_)
  {
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());

    // if geometry is the same, just copy and quit
    if (geom == other.geom_) {
      *this = other;
      return;
    }

    // otherwise, different geometry, do resolution change
    eckit::LocalConfiguration conf;
    conf.set("local interpolator type", "oops unstructured grid interpolator");
    const oops::GeometryData sourceGeom(other.geom_.functionSpace(), other.geom_.fields(),
                                        other.geom_.levelsAreTopDown(), other.geom_.getComm());
    oops::GlobalInterpolator interp(conf, sourceGeom, geom_.functionSpace(), geom.getComm());
    interp.apply(other.fieldSet_, fieldSet_);

    // TODO(Travis) There is a possibility of missing values if the land masks
    // do not match, handle this somehow?

    // TODO(travis) handle a change of resolution in the vertical, someday

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
    soca_state_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    *this = other;
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
    ASSERT(geom_ == rhs.geom_);
    ASSERT(vars_ == rhs.vars_);
    time_ = rhs.time_;
    util::copyFieldSet(rhs.fieldSet_, fieldSet_);
    return *this;
  }

  // -----------------------------------------------------------------------------
  /// Rotations
  // -----------------------------------------------------------------------------
  void State::rotate2north(const oops::Variables & u, const oops::Variables & v) {
    Log::trace() << "State::State rotate from logical to geographical North." << std::endl;

    ASSERT(u.size() == v.size());
    for (size_t n = 0; n < u.size(); n++) {
      const std::string & uName = u[n].name();
      const std::string & vName = v[n].name();
      if (!vars_.has(uName) || !vars_.has(vName)) {
        throw eckit::UserError("State variables " + uName + " or " + vName + " not found.");
      } else {
        Log::info() << "rotating variables " << uName << " and " << vName << std::endl;
      }

      atlas::Field & uField = fieldSet_.field(uName);
      atlas::Field & vField = fieldSet_.field(vName);
      auto uView = atlas::array::make_view<double, 2>(uField);
      auto vView = atlas::array::make_view<double, 2>(vField);
      const auto & ghostView = atlas::array::make_view<int, 1>(uField.functionspace().ghost());
      const auto & cosView = atlas::array::make_view<double, 2>(geom_.fields().field("cos_rot"));
      const auto & sinView = atlas::array::make_view<double, 2>(geom_.fields().field("sin_rot"));

      for (size_t i = 0; i < uField.shape(0); ++i) {
        if (ghostView(i)) continue;

        for (size_t j = 0; j < uField.shape(1); ++j) {
          double uOrig = uView(i, j);
          double vOrig = vView(i, j);
          uView(i, j) = uOrig * cosView(i, 0) + vOrig * sinView(i, 0);
          vView(i, j) = -uOrig * sinView(i, 0) + vOrig * cosView(i, 0);
        }
      }
      uField.set_dirty();
      vField.set_dirty();
    }
  }

  // -----------------------------------------------------------------------------

  void State::rotate2grid(const oops::Variables & u, const oops::Variables & v) {
    Log::trace() << "State::State rotate from geographical to logical North." << std::endl;
    ASSERT(u.size() == v.size());
    for (size_t n = 0; n < u.size(); n++) {
      const std::string & uName = u[n].name();
      const std::string & vName = v[n].name();
      if (!vars_.has(uName) || !vars_.has(vName)) {
        throw eckit::UserError("State variables " + uName + " or " + vName + " not found.");
      } else {
        Log::info() << "rotating variables " << uName << " and " << vName << std::endl;
      }

      atlas::Field & uField = fieldSet_.field(uName);
      atlas::Field & vField = fieldSet_.field(vName);
      auto uView = atlas::array::make_view<double, 2>(uField);
      auto vView = atlas::array::make_view<double, 2>(vField);
      const auto & ghostView = atlas::array::make_view<int, 1>(uField.functionspace().ghost());
      const auto & cosView = atlas::array::make_view<double, 2>(geom_.fields().field("cos_rot"));
      const auto & sinView = atlas::array::make_view<double, 2>(geom_.fields().field("sin_rot"));

      for (size_t i = 0; i < uField.shape(0); ++i) {
        if (ghostView(i)) continue;

        for (size_t j = 0; j < uField.shape(1); ++j) {
          double uOrig = uView(i, j);
          double vOrig = vView(i, j);
          uView(i, j) = uOrig * cosView(i, 0) - vOrig * sinView(i, 0);
          vView(i, j) = uOrig * sinView(i, 0) + vOrig * cosView(i, 0);
        }
      }
      uField.set_dirty();
      vField.set_dirty();
    }
  }

  // -----------------------------------------------------------------------------
  /// Staggered grid interpolation
  // -----------------------------------------------------------------------------
  void State::tohgrid(const oops::Variables & u, const oops::Variables & v) {
    Log::trace() << "State::State interpolate vector to h-grid." << std::endl;
    soca_state_tohgrid_f90(toFortran());
  }

  // -----------------------------------------------------------------------------

  void State::tocgrid(const oops::Variables & u, const oops::Variables & v) {
    Log::trace() << "State::State interpolate vector to c-grid. NOT IMPLEMENTED" << std::endl;
  }

  // -----------------------------------------------------------------------------
  /// Interactions with Increments
  // -----------------------------------------------------------------------------
  State & State::operator+=(const Increment & dx) {
    ASSERT(validTime() == dx.validTime());

    // Interpolate increment to analysis grid only if needed
    std::shared_ptr<const Increment> dx_interp;
    if (geom_ != dx.geometry()) {
      dx_interp = std::make_shared<Increment>(geom_, dx);
    } else {
      dx_interp.reset(&dx, [](const Increment*) {});  // don't delete original dx!
    }

    // Add increment to background state
    // NOTE: if the land masks are not carefully constructed, this can
    // result in MISSING_VALUEs in the increment trying to be added to the state.
    // TODO(travis) issue a warning if this happens? Fix this deeper down in the
    // increment side? In the meantime we can just ignore the missing values

    const auto missing = util::missingValue<double>();
    for (const auto & src : dx_interp->fieldSet()) {
      const auto v_src = atlas::array::make_view<double, 2>(src);
      auto & dst = fieldSet_.field(src.name());
      auto v_dst = atlas::array::make_view<double, 2>(dst);
      for (size_t i = 0; i < src.shape(0); ++i) {
        for (size_t j = 0; j < src.shape(1); ++j) {
          if (v_src(i, j) == missing) continue;
          v_dst(i, j) += v_src(i, j);
        }
      }
      dst.set_dirty(dst.dirty() || src.dirty());
    }

    return *this;
  }

  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------
  void State::read(const eckit::Configuration & files) {
    Log::trace() << "State::State read started." << std::endl;
    util::DateTime * dtp = &time_;
    soca_state_read_file_f90(toFortran(), &files, &dtp);

    fieldSet_.set_dirty();  // just in case, i don't trust the fortran code

    Log::trace() << "State::State read done." << std::endl;
  }

  // -----------------------------------------------------------------------------

  void State::write(const eckit::Configuration & files) const {
    const util::DateTime * dtp = &time_;
    soca_state_write_file_f90(toFortran(), &files, &dtp);
  }

  // -----------------------------------------------------------------------------

  void State::updateFields(const oops::Variables & vars) {
    // remove fields from the fieldset that are no longer in vars
    atlas::FieldSet orig = util::shareFields(fieldSet_);
    fieldSet_.clear();
    for (const auto & v : vars) {
      if (orig.has(v.name())) {
        fieldSet_.add(orig.field(v.name()));
      }
    }

    // update new vars
    vars_ = vars;
    soca_state_update_fields_f90(toFortran(), vars_);
  }

  // -----------------------------------------------------------------------------
  /// Logarithmic and exponential transformations
  // -----------------------------------------------------------------------------

  void State::logtrans(const oops::Variables & trvar) {
    Log::trace() << "State::State apply logarithmic transformation." << std::endl;

    double minVal = 1.0e-6;
    for (const auto & var : trvar) {
      const std::string & varName = var.name();
      if (!vars_.has(varName)) {
        throw eckit::UserError("State variable " + varName + " not found in State.");
      } else {
        Log::info() << "transforming variable "  << varName << std::endl;
      }

      auto & field = fieldSet_.field(varName);
      auto view = atlas::array::make_view<double, 2>(field);
      for (size_t i = 0; i < field.shape(0); ++i) {
        for (size_t j = 0; j < field.shape(1); ++j) {
          view(i, j) = std::log(view(i, j)+minVal);
        }
      }
    }
  }

  // -----------------------------------------------------------------------------

  void State::expontrans(const oops::Variables & trvar) {
    Log::trace() << "State::State apply exponential transformation." << std::endl;

    double minVal = 1.0e-6;
    for (const auto & var : trvar) {
      const std::string & varName = var.name();
      if (!vars_.has(varName)) {
        throw eckit::UserError("State variable " + varName + " not found in State.");
      } else {
        Log::info() << "transforming variable "  << varName << std::endl;
      }

      auto & field = fieldSet_.field(varName);
      auto view = atlas::array::make_view<double, 2>(field);
      for (size_t i = 0; i < field.shape(0); ++i) {
        for (size_t j = 0; j < field.shape(1); ++j) {
          view(i, j) = std::exp(view(i, j))-minVal;
        }
      }
    }
  }

  // -----------------------------------------------------------------------------

}  // namespace soca
