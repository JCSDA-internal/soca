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
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_);
    zero();
    Log::trace() << "Increment constructed." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Geometry & geom, const Increment & other)
    : Fields(geom, other.vars_, other.time_)
  {
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_);
    soca_increment_change_resol_f90(toFortran(), other.keyFlds_);
    Log::trace() << "Increment constructed from other." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other, const bool copy)
    : Fields(other.geom_, other.vars_, other.time_)
  {
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_);
    if (copy) {
      soca_increment_copy_f90(toFortran(), other.toFortran());
    } else {
      zero();
    }
    Log::trace() << "Increment copy-created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other)
    : Fields(other.geom_, other.vars_, other.time_)
  {
    soca_increment_create_f90(keyFlds_, geom_.toFortran(), vars_);
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
    atlas::FieldSet fs1, fs2, fs3;
    x1_at_geomres.toFieldSet(fs2); x2_at_geomres.toFieldSet(fs3);
    fs1 = util::copyFieldSet(fs2);
    util::subtractFieldSets(fs1, fs3);
    fromFieldSet(fs1);
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
    atlas::FieldSet fs1, fs2; toFieldSet(fs1); dx.toFieldSet(fs2);
    util::addFieldSets(fs1, fs2);
    fromFieldSet(fs1);
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator-=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    atlas::FieldSet fs1, fs2; toFieldSet(fs1); dx.toFieldSet(fs2);
    util::subtractFieldSets(fs1, fs2);
    fromFieldSet(fs1);
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator*=(const double & zz) {
    atlas::FieldSet fs1; toFieldSet(fs1);
    util::multiplyFieldSet(fs1, zz);
    fromFieldSet(fs1);
    return *this;
  }
  // -----------------------------------------------------------------------------
  void Increment::ones() {
    atlas::FieldSet fs; toFieldSet(fs);
    for (auto & field : fs) {
      auto view = atlas::array::make_view<double, 2>(field);
      view.assign(1.0);
    }
    fromFieldSet(fs);
  }
  // -----------------------------------------------------------------------------
  void Increment::zero() {
    atlas::FieldSet fs1; toFieldSet(fs1);
    util::zeroFieldSet(fs1);
    fromFieldSet(fs1);
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
  void Increment::axpy(const double & zz, const Increment & dx,
                       const bool check) {
    ASSERT(!check || validTime() == dx.validTime());
    atlas::FieldSet fs1, fs2, fs3; toFieldSet(fs1); dx.toFieldSet(fs2);
    fs3 = util::copyFieldSet(fs2);
    util::multiplyFieldSet(fs3, zz);
    util::addFieldSets(fs1, fs3);
    fromFieldSet(fs1);
  }
  // -----------------------------------------------------------------------------
  void Increment::accumul(const double & zz, const State & xx) {
    atlas::FieldSet fs1, fs2, fs3; toFieldSet(fs1); xx.toFieldSet(fs2);
    fs3 = util::copyFieldSet(fs2);
    util::multiplyFieldSet(fs3, zz);
    util::addFieldSets(fs1, fs3);
    fromFieldSet(fs1);
  }
  // -----------------------------------------------------------------------------
  void Increment::schur_product_with(const Increment & dx) {
    atlas::FieldSet fs1, fs2; toFieldSet(fs1); dx.toFieldSet(fs2);
    util::multiplyFieldSets(fs1, fs2);
    fromFieldSet(fs1);
  }
  // -----------------------------------------------------------------------------
  double Increment::dot_product_with(const Increment & other) const {
    atlas::FieldSet fs1, fs2; toFieldSet(fs1); other.toFieldSet(fs2);
    return util::dotProductFieldSets(fs1, fs2, fs1.field_names(), geom_.getComm());
  }
  // -----------------------------------------------------------------------------
  void Increment::random() {
    soca_increment_random_f90(toFortran());
  }

  // -----------------------------------------------------------------------------
  oops::LocalIncrement Increment::getLocal(const GeometryIterator & iter) const {
    atlas::FieldSet fs; toFieldSet(fs);

    ASSERT(geom_.IteratorDimension() == 2);  // changed need to be made here for 3D
    std::vector<int> varlens(vars_.size());

    // count space needed
    size_t idx = 0;
    for (const auto & var : vars_.variables()) {
      varlens[idx++] = fs.field(var).shape(1);
    }
    size_t totalLen = std::accumulate(varlens.begin(), varlens.end(), 0);

    // fill in vector
    std::vector<double> values;
    values.reserve(totalLen);
    for (const auto & var : vars_.variables()) {
      const auto & view = atlas::array::make_view<double, 2>(fs.field(var));
      for (size_t lvl = 0; lvl < view.shape(1); lvl++) {
        values.push_back(view(iter.i(), lvl));
      }
    }
    ASSERT(values.size() == totalLen);

    return oops::LocalIncrement(vars_, values, varlens);
  }

  // -----------------------------------------------------------------------------
  void Increment::setLocal(const oops::LocalIncrement & values, const GeometryIterator & iter) {
    atlas::FieldSet fs; toFieldSet(fs);
    ASSERT(geom_.IteratorDimension() == 2);  // changes need to be made here for 3D
    const std::vector<double> & vals = values.getVals();
    size_t idx = 0;
    for (const auto & var : vars_.variables()) {
      auto view = atlas::array::make_view<double, 2>(fs.field(var));
      for (size_t lvl = 0; lvl < view.shape(1); lvl++) {
        view(iter.i(), lvl) = vals[idx++];
      }
    }
    ASSERT(idx == vals.size());
    fromFieldSet(fs);
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
    throw eckit::NotImplemented("soca::Increment::rmsByLevel not implemented yet",
                                Here());
  }

  // -----------------------------------------------------------------------------

  void Increment::updateFields(const oops::Variables & vars) {
    // Update local variables
    vars_ = vars;
    // Update field data
    soca_increment_update_fields_f90(toFortran(), vars_);
  }

// -----------------------------------------------------------------------------

  void Increment::toFieldSet(atlas::FieldSet &fs) const {
    soca_increment_to_fieldset_f90(toFortran(), vars_, fs.get());
  }

// -----------------------------------------------------------------------------

  void Increment::fromFieldSet(const atlas::FieldSet &fs) {
    soca_increment_from_fieldset_f90(toFortran(), vars_, fs.get());
  }

// -----------------------------------------------------------------------------

}  // namespace soca
