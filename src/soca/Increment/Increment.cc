/*
 * (C) Copyright 2017-2024 UCAR
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

#include "oops/base/GeometryData.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/generic/GlobalInterpolator.h"
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
    // For now, creation of the fields and their accompanying metadata is done on
    // the Fortran side. This will be moved to the C++ side at a later date.
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_, fieldSet_.get());
    zero();

    Log::trace() << "Increment constructed." << std::endl;
  }

  // -----------------------------------------------------------------------------
  // Resolution change
  Increment::Increment(const Geometry & geom, const Increment & other, const bool ad)
    : Increment(geom, other.vars_, other.time_)
  {
    Log::trace() << "Increment resolution change." << std::endl;

    // same geometry, just copy and quit
    if (geom == other.geom_) {
      soca_increment_copy_f90(toFortran(), other.toFortran());
      return;
    }

    // otherwise, different geometry, do resolution change
    eckit::LocalConfiguration conf;
    conf.set("local interpolator type", "oops unstructured grid interpolator");
    atlas::FieldSet otherFset, selfFset;
    other.toFieldSet(otherFset);
    if (ad) {
      // adjoint interpolation
      const oops::GeometryData sourceGeom(geom_.functionSpace(), geom_.fields(),
                                          geom_.levelsAreTopDown(), geom_.getComm());
      oops::GlobalInterpolator interp(conf, sourceGeom,
                                      other.geom_.functionSpace(), geom.getComm());
      interp.applyAD(selfFset, otherFset);
    } else {
      // interpolation
      const oops::GeometryData sourceGeom(other.geom_.functionSpace(), other.geom_.fields(),
                                          other.geom_.levelsAreTopDown(), other.geom_.getComm());
      oops::GlobalInterpolator interp(conf, sourceGeom, geom_.functionSpace(), geom.getComm());
      interp.apply(otherFset, selfFset);
    }
    fromFieldSet(selfFset);

    // TODO(Travis) There is a possibility of missing values if the land masks
    // do not match, handle this somehow?

    // TODO(travis) handle a change of resolution in the vertical, someday

    Log::trace() << "soca::Increment resolution change DONE." << std::endl;
  }

  // -----------------------------------------------------------------------------

  Increment::Increment(const Increment & other, const bool copy)
    : Increment(other.geom_, other.vars_, other.time_)
  {
    if (copy) {
      soca_increment_copy_f90(toFortran(), other.toFortran());
    }
    Log::trace() << "Increment copy-created." << std::endl;
  }

  // -----------------------------------------------------------------------------

  Increment::Increment(const Increment & other)
    : Increment(other.geom_, other.vars_, other.time_)
  {
    soca_increment_copy_f90(toFortran(), other.toFortran());
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
    util::copyFieldSet(fs1, fieldSet_);
    util::subtractFieldSets(fieldSet_, fs2);
  }

  // -----------------------------------------------------------------------------

  Increment & Increment::operator=(const Increment & rhs) {
    time_ = rhs.time_;
    soca_increment_copy_f90(toFortran(), rhs.toFortran());
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

  void Increment::dirac(const eckit::Configuration & config) {
    soca_increment_dirac_f90(toFortran(), &config);
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
    accumul(zz, dx);
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
    // TODO(travis) use the built-in random. I didn't want to do it with
    // this PR because it would have changed answers.
    soca_increment_random_f90(toFortran());
  }

  // -----------------------------------------------------------------------------

  oops::LocalIncrement Increment::getLocal(const GeometryIterator & iter) const {
    ASSERT(geom_.IteratorDimension() == 2);  // Change will be needed for 3D.
                                             // We don't use 3D right now.
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
    soca_increment_read_file_f90(toFortran(), &files, &dtp);
  }

  // -----------------------------------------------------------------------------

  void Increment::write(const eckit::Configuration & files) const {
    const util::DateTime * dtp = &time_;
    soca_increment_write_file_f90(toFortran(), &files, &dtp);
  }

  // -----------------------------------------------------------------------------

  void Increment::horiz_scales(const eckit::Configuration & config) {
    soca_increment_horiz_scales_f90(toFortran(), &config);
    Log::trace() << "Horiz decorrelation length scales computed." << std::endl;
  }

  // -----------------------------------------------------------------------------

  void Increment::vert_scales(const double & vert) {
    soca_increment_vert_scales_f90(toFortran(), vert);
    Log::trace() << "Vert decorrelation length scales computed." << std::endl;
  }

  // -----------------------------------------------------------------------------

  std::vector<double> Increment::rmsByLevel(const std::string & varname) const {
    throw eckit::NotImplemented("soca::Increment::rmsByLevel not implemented yet", Here());
  }

  // -----------------------------------------------------------------------------

  void Increment::updateFields(const oops::Variables & vars) {
    vars_ = vars;
    soca_increment_update_fields_f90(toFortran(), vars_);
  }

  // -----------------------------------------------------------------------------

}  // namespace soca
