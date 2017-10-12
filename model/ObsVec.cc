
#include <math.h>

#include "util/Logger.h"

#include "model/ObsVec.h"
#include "model/ObsSpace.h"
#include "model/Fortran.h"

namespace mom5cice5 {
// -----------------------------------------------------------------------------
ObsVec::ObsVec(const ObsSpace & obsdb)
  : obsdb_(obsdb), keyOvec_(0)
{
  mom5cice5_obsvec_setup_f90(keyOvec_, obsdb.nout(), obsdb.nobs());
}
// -----------------------------------------------------------------------------
ObsVec::ObsVec(const ObsVec & other, const bool copy)
  : obsdb_(other.obsdb_), keyOvec_(0) {
  mom5cice5_obsvec_clone_f90(other.keyOvec_, keyOvec_);
  if (copy) {
    mom5cice5_obsvec_assign_f90(keyOvec_, other.keyOvec_);
  } else {
    mom5cice5_obsvec_zero_f90(keyOvec_);
  }
}
// -----------------------------------------------------------------------------
ObsVec::~ObsVec() {
  mom5cice5_obsvec_delete_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  mom5cice5_obsvec_assign_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator*= (const double & zz) {
  mom5cice5_obsvec_mul_scal_f90(keyOvec_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator+= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  mom5cice5_obsvec_add_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator-= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  mom5cice5_obsvec_sub_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator*= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  mom5cice5_obsvec_mul_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
ObsVec & ObsVec::operator/= (const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  mom5cice5_obsvec_div_f90(keyOvec_, keyOvecRhs);
  return *this;
}
// -----------------------------------------------------------------------------
void ObsVec::zero() {
  mom5cice5_obsvec_zero_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVec::axpy(const double & zz, const ObsVec & rhs) {
  const int keyOvecRhs = rhs.keyOvec_;
  mom5cice5_obsvec_axpy_f90(keyOvec_, zz, keyOvecRhs);
}
// -----------------------------------------------------------------------------
void ObsVec::invert() {
  mom5cice5_obsvec_invert_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
void ObsVec::random() {
  mom5cice5_obsvec_random_f90(keyOvec_);
}
// -----------------------------------------------------------------------------
double ObsVec::dot_product_with(const ObsVec & other) const {
  const int keyOvecOther = other.keyOvec_;
  double zz;
  mom5cice5_obsvec_dotprod_f90(keyOvec_, keyOvecOther, zz);
  return zz;
}
// -----------------------------------------------------------------------------
double ObsVec::rms() const {
  double zz;
  mom5cice5_obsvec_dotprod_f90(keyOvec_, keyOvec_, zz);
  int iobs;
  mom5cice5_obsvec_nobs_f90(keyOvec_, iobs);
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
  mom5cice5_obsvec_minmaxavg_f90(keyOvec_, zmin, zmax, zavg);
  os << obsdb_.obsname() << " nobs= " << size()
     << " Min=" << zmin << ", Max=" << zmax << ", Average=" << zavg;
}
// -----------------------------------------------------------------------------
unsigned int ObsVec::size() const {
  int iobs;
  mom5cice5_obsvec_nobs_f90(keyOvec_, iobs);
  unsigned int nobs(iobs);
  return nobs;
}
// -----------------------------------------------------------------------------
}  // namespace mom5cice5
