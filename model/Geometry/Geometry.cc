
#include "util/Logger.h"
#include "model/Geometry/Geometry.h"
#include "model/Fortran.h"
#include "eckit/config/Configuration.h"
#include <string>

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration & conf) {
    const eckit::Configuration * configc = &conf;
    std::cout << "***************************** constructor 1" << std::endl;
    soca_geo_setup_f90(keyGeom_, &configc);
  }
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other) {
    const int key_geo = other.keyGeom_;
    std::cout << "***************************** constructor 2 (clone)" << std::endl;    
    soca_geo_clone_f90(key_geo, keyGeom_);
  }
  // -----------------------------------------------------------------------------
  Geometry::~Geometry() {
    soca_geo_delete_f90(keyGeom_);
  }

  // -----------------------------------------------------------------------------          /**/   
  void Geometry::print(std::ostream & os) const {
    int nx;
    int ny;
    int nzo;
    int nzi;
    int ncat;
    soca_geo_info_f90(keyGeom_, nx, ny, nzo, nzi, ncat);
    os << "nx = " << nx << ", ny = " << ny;
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
