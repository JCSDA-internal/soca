

#include "model/Gom.h"

#include "util/Logger.h"
#include "model/ObsSpace.h"
#include "model/Observation.h"
#include "model/LinearObsOp.h"
#include "model/Fortran.h"
#include "model/Variables.h"

namespace mom5cice5 {
  class Geometry;

// -----------------------------------------------------------------------------
Gom::Gom(const ObsSpace & obsdb, const Variables & var,
             const util::DateTime & t1, const util::DateTime & t2,
             const Geometry &) {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  mom5cice5_obsdb_getgom_f90(obsdb.toFortran(), obsdb.obsname().size(), obsdb.obsname().c_str(),
                      var.toFortran(), &p1, &p2, keyGom_);
}
// -----------------------------------------------------------------------------
Gom::~Gom() {
  mom5cice5_gom_delete_f90(keyGom_);
}
// -----------------------------------------------------------------------------
// Gom & Gom::operator= (const Gom & rhs) {
//   const F90ovec * f90rhs = rhs.f90gom_;
//   mom5cice5_gom_assign_f90(keyGom_, &f90rhs);
//   return *this;
// }
// // -----------------------------------------------------------------------------
// Gom & Gom::operator*= (const double & zz) {
//   mom5cice5_gom_mul_scal_f90(keyGom_, zz);
//   return *this;
// }
// // -----------------------------------------------------------------------------
// Gom & Gom::operator+= (const Gom & rhs) {
//   const F90ovec * f90rhs = rhs.f90gom_;
//   mom5cice5_gom_add_f90(keyGom_, &f90rhs);
//   return *this;
// }
// // -----------------------------------------------------------------------------
// Gom & Gom::operator-= (const Gom & rhs) {
//   const F90ovec * f90rhs = rhs.f90gom_;
//   mom5cice5_gom_sub_f90(keyGom_, &f90rhs);
//   return *this;
// }
// // -----------------------------------------------------------------------------
// Gom & Gom::operator*= (const Gom & rhs) {
//   const F90ovec * f90rhs = rhs.f90gom_;
//   mom5cice5_gom_mul_f90(keyGom_, &f90rhs);
//   return *this;
// }
// // -----------------------------------------------------------------------------
// Gom & Gom::operator/= (const Gom & rhs) {
//   const F90ovec * f90rhs = rhs.f90gom_;
//   mom5cice5_gom_div_f90(keyGom_, &f90rhs);
//   return *this;
// }
// -----------------------------------------------------------------------------
void Gom::zero() {
  mom5cice5_gom_zero_f90(keyGom_);
}
// -----------------------------------------------------------------------------
// void Gom::random() {
//   mom5cice5_gom_random_f90(keyGom_);
// }
// -----------------------------------------------------------------------------
double Gom::dot_product_with(const Gom & other) const {
  double zz;
  mom5cice5_gom_dotprod_f90(keyGom_, other.toFortran(), zz);
  return zz;
}
// -----------------------------------------------------------------------------
void Gom::print(std::ostream & os) const {
  int nn;
  double zmin, zmax, zavg;
  mom5cice5_gom_minmaxavg_f90(keyGom_, nn, zmin, zmax, zavg);
  os << " nobs= " << nn << " Min=" << zmin << ", Max=" << zmax << ", RMS=" << zavg;
}
// -----------------------------------------------------------------------------
}  // namespace mom5cice5
