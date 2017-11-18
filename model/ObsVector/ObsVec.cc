
#include <math.h>

#include "util/Logger.h"

#include "model/ObsVector/ObsVec.h"
#include "model/ObsSpace/ObsSpace.h"
#include "model/Fortran.h"

namespace soca {
// -----------------------------------------------------------------------------
ObsVec::ObsVec(const ObsSpace & obsdb)
  : obsdb_(obsdb), keyOvec_(0)
{
  soca_obsvec_setup_f90(keyOvec_, obsdb.nout(), obsdb.nobs());
}
// -----------------------------------------------------------------------------
ObsVec::ObsVec(const ObsVec & other, const bool copy)
  : obsdb_(other.obsdb_), keyOvec_(0) {
  soca_obsvec_clone_f90(other.keyOvec_, keyOvec_);
  if (copy) {
    soca_obsvec_assign_f90(keyOvec_, other.keyOvec_);
  } else {
    soca_obsvec_zero_f90(keyOvec_);
  }
}
// -----------------------------------------------------------------------------
ObsVec::~ObsVec() {
  soca_obsvec_delete_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  soca_obsvec_assign_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator*= (const double & zz) {
  soca_obsvec_mul_scal_f90(keyOvec_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator+= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  soca_obsvec_add_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator-= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  soca_obsvec_sub_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator*= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  soca_obsvec_mul_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator/= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  soca_obsvec_div_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVec::zero() {
  soca_obsvec_zero_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVec::axpy(const double & zz, const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  soca_obsvec_axpy_f90(keyOvec_, zz, keyOvecRhs);
}
// -----------------------------------------------------------------------------
void ObsVec::invert() {
  soca_obsvec_invert_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVec::random() {
  soca_obsvec_random_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
double ObsVec::dot_product_with(const ObsVec & other) const {
  const int keyOvecOther = other.keyOvec_;
  double zz;
  soca_obsvec_dotprod_f90(keyOvec_, keyOvecOther, zz);
  return zz;
}
// -----------------------------------------------------------------------------
double ObsVec::rms() const {
  double zz;
  soca_obsvec_dotprod_f90(keyOvec_, keyOvec_, zz);
  int iobs;
  soca_obsvec_nobs_f90(keyOvec_, iobs);
  zz = sqrt(zz/iobs);
  return zz;
}
// -----------------------------------------------------------------------------
void ObsVec::read(const std::string & name) {
  obsdb_.getdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVec::save(const std::string & name) const {
  obsdb_.putdb(name, keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVec::print(std::ostream & os) const {
  double zmin, zmax, zavg;
  soca_obsvec_minmaxavg_f90(keyOvec_, zmin, zmax, zavg);
  os << obsdb_.obsname() << " nobs= " << size()
     << " Min=" << zmin << ", Max=" << zmax << ", Average=" << zavg;
}
// -----------------------------------------------------------------------------
unsigned int ObsVec::size() const {
  int iobs;
  soca_obsvec_nobs_f90(keyOvec_, iobs);
  unsigned int nobs(iobs);
  return nobs;
}
// -----------------------------------------------------------------------------
}  // namespace soca
