/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/Fields/Fields.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/base/Variables.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"
#include "model/Fortran.h"
#include "model/Geometry/Geometry.h"
#include "model/Variables/Variables.h"

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  Fields::Fields(const Geometry & geom, const oops::Variables & vars,
		 const util::DateTime & time):
    geom_(new Geometry(geom)), vars_(vars), time_(time)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
  }
  // -----------------------------------------------------------------------------
  Fields::Fields(const Fields & other, const bool copy)
    : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
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
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
    soca_field_copy_f90(keyFlds_, other.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  Fields::Fields(const Fields & other, const Geometry & geom)
    : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
    soca_field_change_resol_f90(keyFlds_, other.keyFlds_);
  }
  // -----------------------------------------------------------------------------
  Fields::Fields(const Fields & other, const oops::Variables & vars)
    : geom_(other.geom_), vars_(vars), time_(other.time_)
  {
    soca_field_create_f90(keyFlds_, geom_->toFortran(), vars_.toFortran());
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
    std::cout << "ot of dirac" << std::endl;
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
  void Fields::interpolate(const ioda::Locations & locs, const oops::Variables & vars, ufo::GeoVaLs & gom) const {
    const eckit::Configuration * conf = &vars.toFortran();    
    soca_field_interp_tl_f90(keyFlds_, locs.toFortran(), &conf, gom.toFortran());
   }
  // -----------------------------------------------------------------------------
  void Fields::interpolateTL(const ioda::Locations & locs, const oops::Variables & vars, ufo::GeoVaLs & gom) const {
    const eckit::Configuration * conf = &vars.toFortran();        
    soca_field_interp_tl_f90(keyFlds_, locs.toFortran(), &conf, gom.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Fields::interpolateAD(const ioda::Locations & locs, const oops::Variables & vars, const ufo::GeoVaLs & gom) {
    const eckit::Configuration * conf = &vars.toFortran();            
    soca_field_interp_ad_f90(keyFlds_, locs.toFortran(), &conf, gom.toFortran());
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
  void Fields::convert_to(oops::UnstructuredGrid & ug) const {
    std::cout << "--------- in convert_to " << std::endl;
    soca_field_convert_to_f90(keyFlds_, ug.toFortran());
  }
  // -----------------------------------------------------------------------------
  void Fields::convert_from(const oops::UnstructuredGrid & ug) {
    soca_field_convert_from_f90(keyFlds_, ug.toFortran());
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
    os << std::endl << "  Resolution = " << nx << ", " << ny
       << ", Fields = " << nf;
    std::vector<double> zstat(3*nf);
    soca_field_gpnorm_f90(keyFlds_, nf, zstat[0]);
    for (int jj = 0; jj < nf; ++jj) {
      os << std::endl << "Min=" << zstat[3*jj] << " Max=" << zstat[3*jj+1] << " RMS=" << zstat[3*jj+2];
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
  
    //int nx = -1;
    //int ny = -1;
    //int nf = -1;
    //int nb = -1;
    //soca_field_sizes_f90(keyFlds_, nx, ny, nf, nb);
    std::cout << "================ nf=" << nf << std::endl;
    bool ok = (nf == 11);   //<---- HARD CODED STUFF ... NEED TO CHANGE
    if (nonlinear) ok = ok; // && (nb == 2);
    return ok;
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
