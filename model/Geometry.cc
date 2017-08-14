
#include "util/Logger.h"
#include "model/Geometry.h"
#include "model/Fortran.h"
#include "eckit/config/Configuration.h"
#include <string>

// -----------------------------------------------------------------------------
namespace mom5cice5 {
// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & conf) {
const eckit::Configuration * configc = &conf;
mom5cice5_geo_setup_f90(keyGeom_, &configc);
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) {
const int key_geo = other.keyGeom_;
mom5cice5_geo_clone_f90(key_geo, keyGeom_);
}
// -----------------------------------------------------------------------------
Geometry::~Geometry() {
mom5cice5_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
std::vector<int> Geometry::getDims() const {
std::vector<int> dims(2);
int nzo;
int nzi;
int ncat;
mom5cice5_geo_info_f90(keyGeom_, dims[0], dims[1], nzo, nzi, ncat);  
return dims;
}
// -----------------------------------------------------------------------------                                                                               
std::vector<double> Geometry::getLats() const {
int nx;
int ny;
int nzo;
int nzi;
int ncat;
int geofld;
mom5cice5_geo_info_f90(keyGeom_, nx, ny, nzo, nzi, ncat);
std::vector<double> lats(nx * ny);
  
geofld = 1; 
mom5cice5_geo_getgeofld_f90(keyGeom_, &lats[0], geofld);

return lats;
}
// -----------------------------------------------------------------------------
std::vector<double> Geometry::getLons() const {
int nx;
int ny;
int nzo;
int nzi;
int ncat;
int geofld;
mom5cice5_geo_info_f90(keyGeom_, nx, ny, nzo, nzi, ncat);
std::vector<double> lons(nx * ny);
  
geofld = 0; 
mom5cice5_geo_getgeofld_f90(keyGeom_, &lons[0], geofld);

return lons;
}
// -----------------------------------------------------------------------------
std::vector<double> Geometry::getLevs() const {
int nx;
int ny;
int nzo;
int nzi;
int ncat;
mom5cice5_geo_info_f90(keyGeom_, nx, ny, nzo, nzi, ncat);
std::vector<double> levs(nzi+1);
for (int jj = 0; jj < nzi+1; ++jj) {
levs[jj] = double(jj);
}
return levs;
}
// -----------------------------------------------------------------------------
std::vector<double> Geometry::getArea() const {
int nx;
int ny;
int nzo;
int nzi;
int ncat;
int geofld;
double a;
mom5cice5_geo_info_f90(keyGeom_, nx, ny, nzo, nzi, ncat);

std::vector<double> mask(nx * ny);
geofld=2; // Get the mask
mom5cice5_geo_getgeofld_f90(keyGeom_, &mask[0], geofld);

std::vector<double> cell_area(nx * ny);
geofld=3; // Get the cell area
mom5cice5_geo_getgeofld_f90(keyGeom_, &cell_area[0], geofld);
a=0.0;
for (int jj = 0; jj < nx*ny; ++jj) {
a = a + mask[jj]*cell_area[jj];
}
std::vector<double> area(nzi+1);
for (int jj = 0; jj < nzi+1; ++jj) {
area[jj] = a;
}
return area;
}
// -----------------------------------------------------------------------------
std::vector<int> Geometry::getMask(const int &) const {
int nx;
int ny;
int nzo;
int nzi;
int ncat;
int geofld;
mom5cice5_geo_info_f90(keyGeom_, nx, ny, nzo, nzi, ncat);
std::vector<double> dmask(nx * ny);

geofld=2; // Get the mask
mom5cice5_geo_getgeofld_f90(keyGeom_, &dmask[0], geofld);

std::vector<int> mask(dmask.begin(), dmask.end());
assert(mask.size() == nx*ny);

return mask;
}
// -----------------------------------------------------------------------------             
void Geometry::print(std::ostream & os) const {
int nx;
int ny;
int nzo;
int nzi;
int ncat;
mom5cice5_geo_info_f90(keyGeom_, nx, ny, nzo, nzi, ncat);
os << "nx = " << nx << ", ny = " << ny;
}
// -----------------------------------------------------------------------------
}  // namespace mom5cice5
