/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iomanip>
#include <vector>

#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/Increment/Increment.h"
#include "soca/Increment/IncrementFortran.h"
#include "soca/State/State.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/GridPoint.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
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
  Increment::Increment(const Geometry & geom, const oops::Variables & vars,
                       const util::DateTime & vt)
    : time_(vt), vars_(vars), geom_(new Geometry(geom))
  {
    soca_increment_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_increment_zero_f90(toFortran());
    Log::trace() << "Increment constructed." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Geometry & geom, const Increment & other)
    : time_(other.time_), vars_(other.vars_), geom_(new Geometry(geom))
  {
    soca_increment_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_increment_change_resol_f90(toFortran(), other.keyFlds_);
    Log::trace() << "Increment constructed from other." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other, const bool copy)
    : time_(other.time_), vars_(other.vars_), geom_(new Geometry(*other.geom_))
  {
    soca_increment_create_f90(keyFlds_, geom_->toFortran(), vars_);
    if (copy) {
      soca_increment_copy_f90(toFortran(), other.toFortran());
    } else {
      soca_increment_zero_f90(toFortran());
    }
    Log::trace() << "Increment copy-created." << std::endl;
  }
  // -----------------------------------------------------------------------------
  Increment::Increment(const Increment & other)
    : time_(other.time_), vars_(other.vars_), geom_(new Geometry(*other.geom_))
  {
    soca_increment_create_f90(keyFlds_, geom_->toFortran(), vars_);
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
    Log::debug() << "Increment:diff incr " << *this << std::endl;
    Log::debug() << "Increment:diff x1 " << x1 << std::endl;
    Log::debug() << "Increment:diff x2 " << x2 << std::endl;
    soca_increment_diff_incr_f90(toFortran(), x1.toFortran(),
                             x2.toFortran());
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
    soca_increment_self_add_f90(toFortran(), dx.toFortran());
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator-=(const Increment & dx) {
    ASSERT(this->validTime() == dx.validTime());
    soca_increment_self_sub_f90(toFortran(), dx.toFortran());
    return *this;
  }
  // -----------------------------------------------------------------------------
  Increment & Increment::operator*=(const double & zz) {
    soca_increment_self_mul_f90(toFortran(), zz);
    return *this;
  }
  // -----------------------------------------------------------------------------
  void Increment::zero() {
    soca_increment_zero_f90(toFortran());
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
    soca_increment_axpy_f90(toFortran(), zz, dx.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Increment::accumul(const double & zz, const State & xx) {
    soca_increment_axpy_f90(toFortran(), zz, xx.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Increment::schur_product_with(const Increment & dx) {
    soca_increment_self_schur_f90(toFortran(), dx.toFortran());
  }
  // -----------------------------------------------------------------------------
  double Increment::dot_product_with(const Increment & other) const {
    double zz;
    soca_increment_dot_prod_f90(toFortran(), other.toFortran(), zz);
    return zz;
  }
  // -----------------------------------------------------------------------------
  void Increment::random() {
    soca_increment_random_f90(toFortran());
  }

  // -----------------------------------------------------------------------------
  oops::GridPoint Increment::getPoint(const GeometryIterator & iter) const {
    int nx, ny, nzo, nzi, ncat, nf;
    soca_increment_sizes_f90(toFortran(), nx, ny, nzo, nzi, ncat, nf);

    std::vector<int> varlens(vars_.size());

    // TODO(Travis) remove the hardcoded variable names
    for (int ii = 0; ii < vars_.size(); ii++) {
      if (vars_[ii] == "tocn") varlens[ii]=nzo;
      else if (vars_[ii] == "socn") varlens[ii]=nzo;
      else if (vars_[ii] == "hocn") varlens[ii]=nzo;
      else if (vars_[ii] == "cicen") varlens[ii]=ncat;
      else if (vars_[ii] == "hicen") varlens[ii]=ncat;
      else if (vars_[ii] == "hsnon") varlens[ii]=ncat;
      else
          varlens[ii] = 1;
    }

    int lenvalues = std::accumulate(varlens.begin(), varlens.end(), 0);
    std::vector<double> values(lenvalues);

    soca_increment_getpoint_f90(keyFlds_, iter.toFortran(), values[0],
                            values.size());

    return oops::GridPoint(vars_, values, varlens);
  }

  // -----------------------------------------------------------------------------
  void Increment::setPoint(const oops::GridPoint & values,
                             const GeometryIterator & iter) {
    const std::vector<double> vals = values.getVals();
    soca_increment_setpoint_f90(toFortran(), iter.toFortran(), vals[0],
                            vals.size());
  }

  // -----------------------------------------------------------------------------
  /// Unstructured grid
  // -----------------------------------------------------------------------------
  void Increment::ug_coord(oops::UnstructuredGrid & ug) const {
    soca_increment_ug_coord_f90(toFortran(), ug.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Increment::field_to_ug(oops::UnstructuredGrid & ug,
                              const int & its) const {
    soca_increment_field_to_ug_f90(toFortran(), ug.toFortran(), its);
  }
  // -----------------------------------------------------------------------------
  void Increment::field_from_ug(const oops::UnstructuredGrid & ug,
                                const int & its) {
    soca_increment_field_from_ug_f90(toFortran(), ug.toFortran(), its);
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
  void Increment::print(std::ostream & os) const {
    os << std::endl << "  Valid time: " << validTime();
    int n0, nf;
    soca_increment_sizes_f90(keyFlds_, n0, n0, n0, n0, n0, nf);
    std::vector<double> zstat(3*nf);
    soca_increment_gpnorm_f90(keyFlds_, nf, zstat[0]);
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

  double Increment::norm() const {
    double zz = 0.0;
    soca_increment_rms_f90(toFortran(), zz);
    return zz;
  }

  // -----------------------------------------------------------------------------

  const util::DateTime & Increment::validTime() const {return time_;}

  // -----------------------------------------------------------------------------

  util::DateTime & Increment::validTime() {return time_;}

  // -----------------------------------------------------------------------------

  void Increment::updateTime(const util::Duration & dt) {time_ += dt;}

  // -----------------------------------------------------------------------------

  boost::shared_ptr<const Geometry> Increment::geometry() const {
    return geom_;
  }
  // -----------------------------------------------------------------------------

}  // namespace soca
