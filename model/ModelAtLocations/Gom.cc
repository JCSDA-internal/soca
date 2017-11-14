

#include "model/ModelAtLocations/Gom.h"
#include "util/Logger.h"
#include "model/ObsSpace/ObsSpace.h"
#include "model/Fortran.h"
#include "model/Variables/Variables.h"

namespace mom5cice5 {
  class Geometry;

  // -----------------------------------------------------------------------------
  Gom::Gom(const ObsSpace & obsdb, const Variables & var,
	   const util::DateTime & t1, const util::DateTime & t2) {
    //const Geometry & ) {
    const util::DateTime * p1 = &t1;
    const util::DateTime * p2 = &t2;
    mom5cice5_obsdb_getgom_f90(obsdb.toFortran(), obsdb.obsname().size(), obsdb.obsname().c_str(),
			       var.toFortran(), &p1, &p2, keyGom_);
  }
  // -----------------------------------------------------------------------------
  Gom::Gom(const eckit::Configuration & config) {
    mom5cice5_gom_create_f90(keyGom_);
    const eckit::Configuration * conf = &config;
    std::cout << " gom_read_file " << std::endl;
    mom5cice5_gom_read_file_f90(keyGom_, &conf);
  }
  // -----------------------------------------------------------------------------  
  Gom::~Gom() {
    mom5cice5_gom_delete_f90(keyGom_);
  }
  // -----------------------------------------------------------------------------
  void Gom::zero() {
    mom5cice5_gom_zero_f90(keyGom_);
  }
  // -----------------------------------------------------------------------------
  void Gom::random() {
    mom5cice5_gom_random_f90(keyGom_);
  }
  // -----------------------------------------------------------------------------
  Gom & Gom::operator*=(const double & zz) {
    mom5cice5_gom_mult_f90(keyGom_, zz);
    return *this;
  }
  // -----------------------------------------------------------------------------  
  double Gom::dot_product_with(const Gom & other) const {    
    double zz;
    //mom5cice5_gom_dotprod_f90(keyGom_, other.toFortran(), zz);
    mom5cice5_gom_dotprod_f90(keyGom_, other.keyGom_, zz);    
    return zz;
  }
  // -----------------------------------------------------------------------------
  void Gom::read(const eckit::Configuration & config) {
    const eckit::Configuration * conf = &config;
    mom5cice5_gom_read_file_f90(keyGom_, &conf);
  }
  // -----------------------------------------------------------------------------
  void Gom::write(const eckit::Configuration & config) const {
    const eckit::Configuration * conf = &config;
    std::cout << "---------------- In gom write ---------" << std::endl;
    mom5cice5_gom_write_file_f90(keyGom_, &conf);
  }
  // -----------------------------------------------------------------------------  
  // -----------------------------------------------------------------------------
  void Gom::print(std::ostream & os) const {
    int nn;
    double zmin, zmax, zavg;
    mom5cice5_gom_minmaxavg_f90(keyGom_, nn, zmin, zmax, zavg);
    os << " nobs= " << nn << " Min=" << zmin << ", Max=" << zmax << ", RMS=" << zavg;
  }
  // -----------------------------------------------------------------------------
}  // namespace mom5cice5
