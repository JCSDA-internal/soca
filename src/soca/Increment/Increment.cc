/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iomanip>
#include <numeric>
#include <vector>

#include "atlas/field.h"

#include "eckit/geometry/Point3.h"

#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/Increment/Increment.h"
#include "soca/Increment/IncrementFortran.h"
#include "soca/State/State.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FieldSetHelpers.h"

#include "ufo/GeoVaLs.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  Increment::Increment(const Geometry & geom, const oops::Variables & vars,
                       const util::DateTime & vt)
    : Fields(geom, vars, vt)
  {
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    syncToFieldset();
    zero();

    Log::trace() << "Increment constructed." << std::endl;
  }

  // -----------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom, const Increment & other)
    : Fields(geom, other.vars_, other.time_)
  {
    other.syncFromFieldset();
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    soca_increment_change_resol_f90(toFortran(), other.keyFlds_);
    syncToFieldset();
    Log::trace() << "Increment constructed from other." << std::endl;
  }

  // -----------------------------------------------------------------------------

  Increment::Increment(const Increment & other, const bool copy)
    : Fields(other.geom_, other.vars_, other.time_)
  {
    other.syncFromFieldset();
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    syncToFieldset();
    if (copy) {
      soca_increment_copy_f90(toFortran(), other.toFortran());
      syncToFieldset();
    } else {
      zero();
    }
    Log::trace() << "Increment copy-created." << std::endl;
  }

  // -----------------------------------------------------------------------------

  Increment::Increment(const Increment & other)
    : Fields(other.geom_, other.vars_, other.time_)
  {
    other.syncFromFieldset();
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    soca_increment_copy_f90(toFortran(), other.toFortran());
    syncToFieldset();

    Log::trace() << "Increment copy-created." << std::endl;
  }

  // -----------------------------------------------------------------------------

  Increment::~Increment() {
    soca_increment_delete_f90(toFortran());
    Log::trace() << "Increment destructed" << std::endl;
  }

  // -----------------------------------------------------------------------------
  /// Basic operators
  // -----------------------------------------------------------------------------
  void Increment::diff(const State & x1, const State & x2) {
    ASSERT(this->validTime() == x1.validTime());
    ASSERT(this->validTime() == x2.validTime());

    State x1_at_geomres(geom_, x1);
    State x2_at_geomres(geom_, x2);
    atlas::FieldSet fs1, fs2;
    x1_at_geomres.toFieldSet(fs1); x2_at_geomres.toFieldSet(fs2);
    fieldSet_ = util::copyFieldSet(fs1);
    util::subtractFieldSets(fieldSet_, fs2);
  }

  // -----------------------------------------------------------------------------

  Increment & Increment::operator=(const Increment & rhs) {
    syncFromFieldset();
    time_ = rhs.time_;
    vars_ = rhs.vars_;
    rhs.syncFromFieldset();
    soca_increment_copy_f90(toFortran(), rhs.toFortran());
    syncToFieldset();
    return *this;
  }

  // -----------------------------------------------------------------------------

  Increment & Increment::operator+=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());

    // note, can't use util::addFieldSets because it doesn't handle a variable
    // being in dx but not being in this (not sure why that is happening. is
    // this a bug in soca?)
    for (const auto & addField : dx.fieldSet_) {
      if (!fieldSet_.has(addField.name())) continue;

      atlas::Field field = fieldSet_.field(addField.name());

      auto view = atlas::array::make_view<double, 2>(field);
      const auto addView = atlas::array::make_view<double, 2>(addField);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
        view(jnode, jlevel) += addView(jnode, jlevel);
        }
      }

      // If either term in the sum is out-of-date, then the result will be out-of-date
      field.set_dirty(field.dirty() || addField.dirty());
    }
    return *this;
  }

  // -----------------------------------------------------------------------------

  Increment & Increment::operator-=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    util::subtractFieldSets(fieldSet_, dx.fieldSet_);
    return *this;
  }

  // -----------------------------------------------------------------------------

  Increment & Increment::operator*=(const double & zz) {
    util::multiplyFieldSet(fieldSet_, zz);
    return *this;
  }

  // -----------------------------------------------------------------------------

  void Increment::ones() {
    for (auto & field : fieldSet_) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(1.0);
    }
  }

  // -----------------------------------------------------------------------------

  void Increment::zero() {
    util::zeroFieldSet(fieldSet_);
  }

  // -----------------------------------------------------------------------------

  void Increment::dirac(const eckit::Configuration & config) {
    syncFromFieldset();
    soca_increment_dirac_f90(toFortran(), &config);
    syncToFieldset();
    Log::trace() << "Increment dirac initialized" << std::endl;
  }

  // -----------------------------------------------------------------------------
  void Increment::zero(const util::DateTime & vt) {
    zero();
    time_ = vt;
  }

  // -----------------------------------------------------------------------------

  void Increment::axpy(const double & zz, const Increment & dx, const bool check) {
    ASSERT(!check || validTime() == dx.validTime());
    atlas::FieldSet fs1;
    fs1 = util::copyFieldSet(dx.fieldSet_);
    util::multiplyFieldSet(fs1, zz);
    util::addFieldSets(fieldSet_, fs1);
  }

  // -----------------------------------------------------------------------------

  void Increment::accumul(const double & zz, const State & xx) {
    atlas::FieldSet fs1, fs2; xx.toFieldSet(fs2);
    fs1 = util::copyFieldSet(fs2);
    util::multiplyFieldSet(fs1, zz);
    util::addFieldSets(fieldSet_, fs1);
  }

  // -----------------------------------------------------------------------------

  void Increment::schur_product_with(const Increment & dx) {
    // note, can't use util::multiplyFieldSets because it doesn't handle a variable
    // being in dx but not being in this (not sure why that is happening. is
    // this a bug in soca?)
    for (const auto & mulField : dx.fieldSet_) {
      if (!fieldSet_.has(mulField.name())) continue;
      // Get field with the same name
      atlas::Field field = fieldSet_.field(mulField.name());
      auto view = atlas::array::make_view<double, 2>(field);
      const auto mulView = atlas::array::make_view<double, 2>(mulField);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= mulView(jnode, jlevel);
        }
      }
      // If either term in the product is out-of-date, then the result will be out-of-date
      field.set_dirty(field.dirty() || mulField.dirty());
    }
  }

  // -----------------------------------------------------------------------------

  double Increment::dot_product_with(const Increment & other) const {
    return util::dotProductFieldSets(fieldSet_, other.fieldSet_,
      fieldSet_.field_names(), geom_.getComm());
  }

  // -----------------------------------------------------------------------------

  void Increment::random() {
    syncFromFieldset();
    soca_increment_random_f90(toFortran());
    syncToFieldset();
  }

  // -----------------------------------------------------------------------------

  oops::LocalIncrement Increment::getLocal(const GeometryIterator & iter) const {
    ASSERT(geom_.IteratorDimension() == 2);  // changed need to be made here for 3D
    std::vector<int> varlens(vars_.size());

    // count space needed
    size_t idx = 0;
    for (const auto & var : vars_.variables()) {
      varlens[idx++] = fieldSet_.field(var).shape(1);
    }
    size_t totalLen = std::accumulate(varlens.begin(), varlens.end(), 0);

    // fill in vector
    std::vector<double> values;
    values.reserve(totalLen);
    for (const auto & var : vars_.variables()) {
      const auto & view = atlas::array::make_view<double, 2>(fieldSet_.field(var));
      for (size_t lvl = 0; lvl < view.shape(1); lvl++) {
        values.push_back(view(iter.i(), lvl));
      }
    }
    ASSERT(values.size() == totalLen);
    return oops::LocalIncrement(vars_, values, varlens);
  }

  // -----------------------------------------------------------------------------

  void Increment::setLocal(const oops::LocalIncrement & values, const GeometryIterator & iter) {
    ASSERT(geom_.IteratorDimension() == 2);  // changes need to be made here for 3D
    const std::vector<double> & vals = values.getVals();
    size_t idx = 0;
    for (const auto & var : vars_.variables()) {
      auto view = atlas::array::make_view<double, 2>(fieldSet_.field(var));
      for (size_t lvl = 0; lvl < view.shape(1); lvl++) {
        view(iter.i(), lvl) = vals[idx++];
      }
    }
    ASSERT(idx == vals.size());
  }

  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------

  void Increment::read(const eckit::Configuration & files) {
    util::DateTime * dtp = &time_;
    syncFromFieldset();
    soca_increment_read_file_f90(toFortran(), &files, &dtp);
    syncToFieldset();
  }

  // -----------------------------------------------------------------------------

  void Increment::write(const eckit::Configuration & files) const {
    const util::DateTime * dtp = &time_;
    syncFromFieldset();
    soca_increment_write_file_f90(toFortran(), &files, &dtp);
  }

  // -----------------------------------------------------------------------------

  void Increment::horiz_scales(const eckit::Configuration & config) {
    syncFromFieldset();
    soca_increment_horiz_scales_f90(toFortran(), &config);
    syncToFieldset();
    Log::trace() << "Horiz decorrelation length scales computed." << std::endl;
  }

  // -----------------------------------------------------------------------------

  void Increment::vert_scales(const double & vert) {
    syncFromFieldset();
    soca_increment_vert_scales_f90(toFortran(), vert);
    syncToFieldset();
    Log::trace() << "Vert decorrelation length scales computed." << std::endl;
  }

  // -----------------------------------------------------------------------------

  std::vector<double> Increment::rmsByLevel(const std::string & varname) const {
    throw eckit::NotImplemented("soca::Increment::rmsByLevel not implemented yet", Here());
  }

  // -----------------------------------------------------------------------------

  void Increment::updateFields(const oops::Variables & vars) {
    syncFromFieldset();
    vars_ = vars;
    soca_increment_update_fields_f90(toFortran(), vars_);
    syncToFieldset();
  }

  // -----------------------------------------------------------------------------

  void Increment::toFieldSet(atlas::FieldSet &fset) const {
    util::copyFieldSet(fieldSet_, fset);
  }

  // -----------------------------------------------------------------------------

  void Increment::fromFieldSet(const atlas::FieldSet &fs) {
    soca_increment_from_fieldset_f90(toFortran(), vars_, fs.get());
    // due to a bug in the metadata being set, i think
    soca_increment_to_fieldset_f90(toFortran(), vars_, fieldSet_.get());
  }

  // -----------------------------------------------------------------------------
  void Increment::syncFromFieldset() const {
    soca_increment_from_fieldset_f90(toFortran(), vars_, fieldSet_.get());
    // due to a bug in the metadata being set, i think
    soca_increment_to_fieldset_f90(toFortran(), vars_, fieldSet_.get());
  }

  // -----------------------------------------------------------------------------

  void Increment::syncToFieldset() const {
    soca_increment_to_fieldset_f90(toFortran(), vars_, fieldSet_.get());
  }

}  // namespace soca
