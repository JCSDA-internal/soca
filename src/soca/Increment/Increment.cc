/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iomanip>
#include <vector>

#include "atlas/field.h"

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
    State x1_at_geomres(*geom_, x1);
    State x2_at_geomres(*geom_, x2);
    soca_increment_diff_incr_f90(toFortran(), x1_at_geomres.toFortran(),
                                              x2_at_geomres.toFortran());
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
  void Increment::ones() {
    soca_increment_ones_f90(toFortran());
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
    soca_increment_accumul_f90(toFortran(), zz, xx.toFortran());
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
  oops::LocalIncrement Increment::getLocal(
                        const GeometryIterator & iter) const {
    int nx, ny, nzo, nzi, ncat, nf;
    soca_increment_sizes_f90(toFortran(), nx, ny, nzo, nzi, ncat, nf);

    std::vector<int> varlens(vars_.size());

    // TODO(Travis) remove the hardcoded variable names
    for (int ii = 0; ii < vars_.size(); ii++) {
      if (vars_[ii] == "tocn") varlens[ii]=nzo;
      else if (vars_[ii] == "socn") varlens[ii]=nzo;
      else if (vars_[ii] == "hocn") varlens[ii]=nzo;
      else if (vars_[ii] == "uocn") varlens[ii]=nzo;
      else if (vars_[ii] == "vocn") varlens[ii]=nzo;
      else if (vars_[ii] == "ssh")  varlens[ii]=1;
      else if (vars_[ii] == "cicen") varlens[ii]=ncat;
      else if (vars_[ii] == "hicen") varlens[ii]=ncat;
      else if (vars_[ii] == "hsnon") varlens[ii]=ncat;
      else if (vars_[ii] == "chl") varlens[ii]=nzo;
      else
          varlens[ii] = 0;
    }

    int lenvalues = std::accumulate(varlens.begin(), varlens.end(), 0);
    std::vector<double> values(lenvalues);

    soca_increment_getpoint_f90(keyFlds_, iter.toFortran(), values[0],
                            values.size());

    return oops::LocalIncrement(vars_, values, varlens);
  }

  // -----------------------------------------------------------------------------
  void Increment::setLocal(const oops::LocalIncrement & values,
                             const GeometryIterator & iter) {
    const std::vector<double> vals = values.getVals();
    soca_increment_setpoint_f90(toFortran(), iter.toFortran(), vals[0],
                            vals.size());
  }
  // -----------------------------------------------------------------------------
  /// ATLAS
  // -----------------------------------------------------------------------------
  void Increment::setAtlas(atlas::FieldSet * afieldset) const {
    soca_increment_set_atlas_f90(toFortran(), geom_->toFortran(), vars_,
                                 afieldset->get());
  }
  // -----------------------------------------------------------------------------
  void Increment::toAtlas(atlas::FieldSet * afieldset) const {
    soca_increment_to_atlas_f90(toFortran(), geom_->toFortran(), vars_,
                                afieldset->get());
  }
  // -----------------------------------------------------------------------------
  void Increment::fromAtlas(atlas::FieldSet * afieldset) {
    soca_increment_from_atlas_f90(toFortran(), geom_->toFortran(), vars_,
                                  afieldset->get());
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
  /// Serialization
  // -----------------------------------------------------------------------------
  size_t Increment::serialSize() const {
    // Field
    size_t nn;
    soca_increment_serial_size_f90(toFortran(), geom_->toFortran(), nn);

    // Magic factor
    nn += 1;

    // Date and time
    nn += time_.serialSize();
    return nn;
  }
  // -----------------------------------------------------------------------------
  constexpr double SerializeCheckValue = -54321.98765;
  void Increment::serialize(std::vector<double> & vect) const {
    // Serialize the field
    size_t nn;
    soca_increment_serial_size_f90(toFortran(), geom_->toFortran(), nn);
    std::vector<double> vect_field(nn, 0);
    vect.reserve(vect.size() + nn + 1 + time_.serialSize());
    soca_increment_serialize_f90(toFortran(), geom_->toFortran(), nn,
                                 vect_field.data());
    vect.insert(vect.end(), vect_field.begin(), vect_field.end());

    // Magic value placed in serialization; used to validate deserialization
    vect.push_back(SerializeCheckValue);

    // Serialize the date and time
    time_.serialize(vect);
  }
  // -----------------------------------------------------------------------------
  void Increment::deserialize(const std::vector<double> & vect,
                              size_t & index) {
    // Deserialize the field

    soca_increment_deserialize_f90(toFortran(), geom_->toFortran(), vect.size(),
                                   vect.data(), index);

    // Use magic value to validate deserialization
    ASSERT(vect.at(index) == SerializeCheckValue);
    ++index;

    // Deserialize the date and time
    time_.deserialize(vect, index);
  }
  // -----------------------------------------------------------------------------

  std::shared_ptr<const Geometry> Increment::geometry() const {
    return geom_;
  }
  // -----------------------------------------------------------------------------

}  // namespace soca
