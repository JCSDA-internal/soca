/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "soca/Fields/Fields.h"
#include "soca/Fields/FieldsFortran.h"
#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/GetValuesTraj/GetValuesTraj.h"

#include "eckit/config/Configuration.h"
#include "oops/base/GridPoint.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  Fields::Fields(const Geometry & geom,
                 const oops::Variables & vars,
                 const util::DateTime & time):
    geom_(new Geometry(geom)), vars_(vars), time_(time)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
  }
  // -----------------------------------------------------------------------------
  Fields::Fields(const Fields & other, const bool copy)
    : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    if (copy) {
      soca_field_copy_f90(keyFlds_, other.keyFlds_);
    } else {
      soca_field_zero_f90(keyFlds_);
    }
  }
  // -----------------------------------------------------------------------------
  Fields::Fields(const Fields & other)
    : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_field_copy_f90(keyFlds_, other.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  Fields::Fields(const Fields & other, const Geometry & geom)
    : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_field_change_resol_f90(keyFlds_, other.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  Fields::Fields(const Fields & other, const oops::Variables & vars)
    : geom_(other.geom_), vars_(vars), time_(other.time_)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_);
    soca_field_copy_f90(keyFlds_, other.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  Fields::~Fields() {
    soca_field_delete_f90(keyFlds_);
  }
  // -----------------------------------------------------------------------------
  Fields & Fields::operator=(const Fields & rhs) {
    soca_field_copy_f90(keyFlds_, rhs.keyFlds_);
    time_ = rhs.time_;
    return *this;
  }
  // -----------------------------------------------------------------------------
  Fields & Fields::operator+=(const Fields & rhs) {
    soca_field_self_add_f90(keyFlds_, rhs.keyFlds_);
    return *this;
  }
  // -----------------------------------------------------------------------------
  Fields & Fields::operator-=(const Fields & rhs) {
    soca_field_self_sub_f90(keyFlds_, rhs.keyFlds_);
    return *this;
  }
  // -----------------------------------------------------------------------------
  Fields & Fields::operator*=(const double & zz) {
    soca_field_self_mul_f90(keyFlds_, zz);
    return *this;
  }
  // -----------------------------------------------------------------------------
  void Fields::zero() {
    soca_field_zero_f90(keyFlds_);
  }
  // -----------------------------------------------------------------------------
  void Fields::zero(const util::DateTime & time) {
    soca_field_zero_f90(keyFlds_);
    time_ = time;
  }
  // -----------------------------------------------------------------------------
  void Fields::dirac(const eckit::Configuration & config) {
    const eckit::Configuration * conf = &config;
    soca_field_dirac_f90(keyFlds_, &conf);
  }
  // -----------------------------------------------------------------------------
  void Fields::axpy(const double & zz, const Fields & rhs) {
    soca_field_axpy_f90(keyFlds_, zz, rhs.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  double Fields::dot_product_with(const Fields & fld2) const {
    double zz;
    soca_field_dot_prod_f90(keyFlds_, fld2.keyFlds_, zz);
    return zz;
  }
  // -----------------------------------------------------------------------------
  void Fields::schur_product_with(const Fields & dx) {
    soca_field_self_schur_f90(keyFlds_, dx.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  void Fields::random() {
    soca_field_random_f90(keyFlds_);
  }
  // -----------------------------------------------------------------------------
  void Fields::getValues(const ufo::Locations & locs,
                         const oops::Variables & vars,
                         ufo::GeoVaLs & gom) const {
    soca_field_interp_nl_f90(keyFlds_, locs.toFortran(), vars,
                             gom.toFortran());
  }

  // -----------------------------------------------------------------------------
  void Fields::getValues(const ufo::Locations & locs,
                         const oops::Variables & vars,
                         ufo::GeoVaLs & gom,
                         const GetValuesTraj & traj) const {
    soca_field_interp_nl_traj_f90(keyFlds_, locs.toFortran(), vars,
                                  gom.toFortran(), traj.toFortran());
  }

  // -----------------------------------------------------------------------------
  void Fields::getValuesTL(const ufo::Locations & locs,
                           const oops::Variables & vars,
                           ufo::GeoVaLs & gom,
                           const GetValuesTraj & traj) const {
    soca_field_interp_tl_f90(keyFlds_, locs.toFortran(), vars,
                             gom.toFortran(), traj.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Fields::getValuesAD(const ufo::Locations & locs,
                           const oops::Variables & vars,
                           const ufo::GeoVaLs & gom,
                           const GetValuesTraj & traj) {
    soca_field_interp_ad_f90(keyFlds_, locs.toFortran(), vars,
                             gom.toFortran(), traj.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Fields::rotate2north() const {
    // rotate2north:
    //     Rotate horizontal vector (u_t,v_t) from logical to true North
    // u_t =   cos(theta) u_l + sin(theta) v_l
    // v_t = - sin(theta) u_l + cos(theta) v_l
    // Assumes u_t and v_t are colocated
    soca_field_rotate2north_f90(keyFlds_);
  }
  // -----------------------------------------------------------------------------
  void Fields::rotate2grid() const {
    // Inverse of rotate2north:
    //     Rotate horizontal vector (u_t,v_t) from true zonal/meridional
    //                           to (u_l,v_l) in local coordinates
    // u_l = cos(theta) u_t - sin(theta) v_t
    // v_l = sin(theta) u_t + cos(theta) v_t
    // Assumes u_t and v_t are colocated
    soca_field_rotate2grid_f90(keyFlds_);
  }
  // -----------------------------------------------------------------------------
  oops::GridPoint Fields::getPoint(const soca::GeometryIterator & iter) {
    int nx, ny, nzo, nzi, ncat, nf;
    soca_field_sizes_f90(toFortran(), nx, ny, nzo, nzi, ncat, nf);

    oops::Variables fieldNames = variables();
    std::vector<int> varlens(fieldNames.size());

    // TODO(Travis) remove the hardcoded variable names
    for (int ii = 0; ii < fieldNames.size(); ii++) {
      if (fieldNames[ii] == "tocn") varlens[ii]=nzo;
      else if (fieldNames[ii] == "socn") varlens[ii]=nzo;
      else if (fieldNames[ii] == "hocn") varlens[ii]=nzo;
      else if (fieldNames[ii] == "cicen") varlens[ii]=ncat;
      else if (fieldNames[ii] == "hicen") varlens[ii]=ncat;
      else if (fieldNames[ii] == "hsnon") varlens[ii]=ncat;
      else
          varlens[ii] = 1;
    }

    int lenvalues = std::accumulate(varlens.begin(), varlens.end(), 0);
    std::vector<double> values(lenvalues);

    soca_field_getpoint_f90(keyFlds_, iter.toFortran(), values[0],
                            values.size());

    return oops::GridPoint(oops::Variables(fieldNames), values, varlens);
  }
  // -----------------------------------------------------------------------------
  void Fields::setPoint(const soca::GeometryIterator & iter,
                const oops::GridPoint & values) {
    const std::vector<double> vals = values.getVals();
    soca_field_setpoint_f90(keyFlds_, iter.toFortran(), vals[0], vals.size());
  }
  // -----------------------------------------------------------------------------
  void Fields::changeResolution(const Fields & other) {
    soca_field_change_resol_f90(keyFlds_, other.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  void Fields::add(const Fields & rhs) {
    soca_field_add_incr_f90(keyFlds_, rhs.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  void Fields::diff(const Fields & x1, const Fields & x2) {
    soca_field_diff_incr_f90(keyFlds_, x1.keyFlds_, x2.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  void Fields::ug_coord(oops::UnstructuredGrid & ug) const {
    soca_field_ug_coord_f90(keyFlds_, ug.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Fields::field_to_ug(oops::UnstructuredGrid & ug,
                           const int & its) const {
    soca_field_field_to_ug_f90(keyFlds_, ug.toFortran(), its);
  }
  // -----------------------------------------------------------------------------
  void Fields::field_from_ug(const oops::UnstructuredGrid & ug,
                             const int & its) {
    soca_field_field_from_ug_f90(keyFlds_, ug.toFortran(), its);
  }
  // -----------------------------------------------------------------------------
  void Fields::read(const eckit::Configuration & config) {
    const eckit::Configuration * conf = &config;
    util::DateTime * dtp = &time_;
    soca_field_read_file_f90(keyFlds_, &conf, &dtp);
  }
  // -----------------------------------------------------------------------------
  void Fields::write(const eckit::Configuration & config) const {
    const eckit::Configuration * conf = &config;
    const util::DateTime * dtp = &time_;
    soca_field_write_file_f90(keyFlds_, &conf, &dtp);
  }
  // -----------------------------------------------------------------------------
  double Fields::norm() const {
    double zz = 0.0;
    soca_field_rms_f90(keyFlds_, zz);
    return zz;
  }
  // -----------------------------------------------------------------------------
  void Fields::print(std::ostream & os) const {
    int nx = -1;
    int ny = -1;
    int nzo = -1;
    int nzi = -1;
    int ncat = -1;
    int nf = -1;
    soca_field_sizes_f90(keyFlds_, nx, ny, nzo, nzi, ncat, nf);
    std::vector<double> zstat(3*nf);
    soca_field_gpnorm_f90(keyFlds_, nf, zstat[0]);
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
  bool Fields::isForModel(bool nonlinear) const {
    int nx = -1;
    int ny = -1;
    int nzo = -1;
    int nzi = -1;
    int ncat = -1;
    int nf = -1;
    soca_field_sizes_f90(keyFlds_, nx, ny, nzo, nzi, ncat, nf);
    bool ok = (nf >= 3);    // <---- HARD CODED STUFF ... NEED TO CHANGE
    if (nonlinear) ok = ok;  // && (nb == 2);
    return ok;
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
