/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Increment/Increment.h"

#include <algorithm>
#include <string>
#include <vector>
#include <utility>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "src/ModelBiasIncrement.h"
#include "src/Covariance/ErrorCovariance.h"
#include "src/Fields/Fields.h"
#include "src/Geometry/Geometry.h"
#include "src/State/State.h"
#include "src/GetValuesTraj/GetValuesTraj.h"
#include "oops/util/Duration.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

using oops::Log;

namespace soca {

  // -----------------------------------------------------------------------------
  /// Constructor, destructor
  // -----------------------------------------------------------------------------
  Increment::Increment(const Geometry & resol, const oops::Variables & vars,
                       const util::DateTime & vt)
    : fields_(new Fields(resol, vars, vt)), stash_()
  {
    fields_->zero();
    Log::trace() << "Increment constructed." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Geometry & resol, const Increment & other)
    : fields_(new Fields(*other.fields_, resol)), stash_()
  {
    Log::trace() << "Increment constructed from other." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other, const bool copy)
    : fields_(new Fields(*other.fields_, copy)), stash_()
  {
    Log::trace() << "Increment copy-created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other)
    : fields_(new Fields(*other.fields_)), stash_()
  {
    Log::trace() << "Increment copy-created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::~Increment() {
    Log::trace() << "Increment destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  /// Basic operators
  // -----------------------------------------------------------------------------
  void Increment::diff(const State & x1, const State & x2) {
    ASSERT(this->validTime() == x1.validTime());
    ASSERT(this->validTime() == x2.validTime());
    Log::debug() << "Increment:diff incr " << *fields_ << std::endl;
    Log::debug() << "Increment:diff x1 " << x1.fields() << std::endl;
    Log::debug() << "Increment:diff x2 " << x2.fields() << std::endl;
    fields_->diff(x1.fields(), x2.fields());
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator=(const Increment & rhs) {
    *fields_ = *rhs.fields_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator+=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    *fields_ += *dx.fields_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator-=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    *fields_ -= *dx.fields_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator*=(const double & zz) {
    *fields_ *= zz;
    return *this;
  }
  // -----------------------------------------------------------------------------
  void Increment::zero() {
    fields_->zero();
  }
  // -----------------------------------------------------------------------------
  void Increment::dirac(const eckit::Configuration & config) {
    fields_->dirac(config);
    Log::trace() << "Increment dirac initialized" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Increment::zero(const util::DateTime & vt) {
    fields_->zero(vt);
  }
  // -----------------------------------------------------------------------------
  void Increment::axpy(const double & zz, const Increment & dx,
                       const bool check) {
    ASSERT(!check || this->validTime() == dx.validTime());
    fields_->axpy(zz, *dx.fields_);
  }
  // -----------------------------------------------------------------------------
  void Increment::accumul(const double & zz, const State & xx) {
    fields_->axpy(zz, xx.fields());
  }
  // -----------------------------------------------------------------------------
  void Increment::schur_product_with(const Increment & dx) {
    fields_->schur_product_with(*dx.fields_);
  }
  // -----------------------------------------------------------------------------
  double Increment::dot_product_with(const Increment & other) const {
    return dot_product(*fields_, *other.fields_);
  }
  // -----------------------------------------------------------------------------
  void Increment::random() {
    fields_->random();
  }
  /// Interpolate to observation location
  // -----------------------------------------------------------------------------
  void Increment::getValuesTL(const ufo::Locations & locs,
                              const oops::Variables & vars,
                              ufo::GeoVaLs & cols,
                              const GetValuesTraj & traj) const {
    Log::debug() << "Increment::interpolateTL fields in" <<
                    *fields_ << std::endl;
    fields_->getValuesTL(locs, vars, cols, traj);
    Log::debug() << "Increment::interpolateTL " << cols << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Increment::getValuesAD(const ufo::Locations & locs,
                              const oops::Variables & vars,
                              const ufo::GeoVaLs & cols,
                              const GetValuesTraj & traj) {
    Log::debug() << "Increment::interpolateAD gom " << cols << std::endl;
    Log::debug() << "Increment::interpolateAD fields in" <<
                    *fields_ << std::endl;
    fields_->getValuesAD(locs, vars, cols, traj);
  }
  // -----------------------------------------------------------------------------
  /// Unstructured grid
  // -----------------------------------------------------------------------------
  void Increment::ug_coord(oops::UnstructuredGrid & ug) const {
    fields_->ug_coord(ug);
  }
  // -----------------------------------------------------------------------------
  void Increment::field_to_ug(oops::UnstructuredGrid & ug,
                              const int & its) const {
    fields_->field_to_ug(ug, its);
  }
  // -----------------------------------------------------------------------------
  void Increment::field_from_ug(const oops::UnstructuredGrid & ug,
                                const int & its) {
    fields_->field_from_ug(ug, its);
  }
  // -----------------------------------------------------------------------------
  /// I/O and diagnostics
  // -----------------------------------------------------------------------------
  void Increment::read(const eckit::Configuration & files) {
    fields_->read(files);
  }
  // -----------------------------------------------------------------------------
  void Increment::write(const eckit::Configuration & files) const {
    fields_->write(files);
  }
  // -----------------------------------------------------------------------------
  void Increment::print(std::ostream & os) const {
    // os << std::endl << "  Valid time: " << validTime();
    // os << *fields_;
  }
  // -----------------------------------------------------------------------------

}  // namespace soca
