

#include "model/ModelAtLocations/Gom.h"
#include "util/Logger.h"
#include "model/ObsSpace/ObsSpace.h"
#include "model/Fortran.h"
#include "model/Variables/Variables.h"

namespace soca {
  class Geometry;

  // -----------------------------------------------------------------------------
  Gom::Gom(const ObsSpace & obsdb, const Variables & var,
	   const util::DateTime & t1, const util::DateTime & t2) {
    //const Geometry & ) {
    const util::DateTime * p1 = &t1;
    const util::DateTime * p2 = &t2;
    soca_obsdb_getgom_f90(obsdb.toFortran(), obsdb.obsname().size(), obsdb.obsname().c_str(),
			       var.toFortran(), &p1, &p2, keyGom_);
  }
  // -----------------------------------------------------------------------------
  Gom::Gom(const eckit::Configuration & config) {
    soca_gom_create_f90(keyGom_);
    const eckit::Configuration * conf = &config;
    std::cout << " gom_read_file " << std::endl;
    soca_gom_read_file_f90(keyGom_, &conf);
  }
  // -----------------------------------------------------------------------------  
  Gom::~Gom() {
    soca_gom_delete_f90(keyGom_);
  }
  // -----------------------------------------------------------------------------
  void Gom::zero() {
    soca_gom_zero_f90(keyGom_);
  }
  // -----------------------------------------------------------------------------
  void Gom::random() {
    soca_gom_random_f90(keyGom_);
  }
  // -----------------------------------------------------------------------------
  Gom & Gom::operator*=(const double & zz) {
    soca_gom_mult_f90(keyGom_, zz);
    return *this;
  }
  // -----------------------------------------------------------------------------  
  double Gom::dot_product_with(const Gom & other) const {    
    double zz;
    //soca_gom_dotprod_f90(keyGom_, other.toFortran(), zz);
    soca_gom_dotprod_f90(keyGom_, other.keyGom_, zz);    
    return zz;
  }
  // -----------------------------------------------------------------------------
  void Gom::read(const eckit::Configuration & config) {
    const eckit::Configuration * conf = &config;
    soca_gom_read_file_f90(keyGom_, &conf);
  }
  // -----------------------------------------------------------------------------
  void Gom::write(const eckit::Configuration & config) const {
    const eckit::Configuration * conf = &config;
    std::cout << "---------------- In gom write ---------" << std::endl;
    soca_gom_write_file_f90(keyGom_, &conf);
  }
  // -----------------------------------------------------------------------------  
  // -----------------------------------------------------------------------------
  void Gom::print(std::ostream & os) const {
    int nn;
    double zmin, zmax, zavg;
    soca_gom_minmaxavg_f90(keyGom_, nn, zmin, zmax, zavg);
    os << " nobs= " << nn << " Min=" << zmin << ", Max=" << zmax << ", RMS=" << zavg;
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
