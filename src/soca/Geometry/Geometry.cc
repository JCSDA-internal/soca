/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/DateTime.h"
#include "soca/Fortran.h"
#include "soca/Geometry/GeometryFortran.h"
#include "soca/Geometry/Geometry.h"
#include "eckit/config/Configuration.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm)
  {
    const eckit::Configuration * configc = &conf;

    // Check if a basic atm  geometry is requiered
    bool init_atm_geom = configc->getBool("atm_geom.init", false);

    // If init = true, create the obs space atmospheric geometry
    if (init_atm_geom)
      {
        // Get Time Bounds
        const util::DateTime bgn = util::DateTime(conf.getString("atm_geom.date_begin"));
        const util::DateTime end = util::DateTime(conf.getString("atm_geom.date_end"));

        // Create the Atmospheric Geometry in Observation Space
        // Weird, but look at it as a time dependent unstructured grid for the atmosphere
        const eckit::LocalConfiguration confatmgeom(conf, "atm_geom.ObsSpace");
        ioda::ObsSpace atmobs( confatmgeom, comm, bgn, end);

        // Get grid size
        int nlocs = atmobs.nlocs();

        // Get locations
        std::vector<double> lats(nlocs);
        std::vector<double> lons(nlocs);
        atmobs.get_db("MetaData", "longitude", nlocs, lons.data());
        atmobs.get_db("MetaData", "latitude", nlocs, lats.data());
        soca_geo_setup_f90(keyGeom_, &configc, nlocs, &lats[0], &lons[0]);
      } else
      {
        int nlocs = 1;
        std::vector<double> lats(nlocs); lats[0] = 0.0;
        std::vector<double> lons(nlocs); lons[0] = 0.0;
        soca_geo_setup_f90(keyGeom_, &configc, nlocs, &lats[0], &lons[0]);
      }
  }
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_) {
    const int key_geo = other.keyGeom_;
    soca_geo_clone_f90(key_geo, keyGeom_);
  }
  // -----------------------------------------------------------------------------
  Geometry::~Geometry() {
    soca_geo_delete_f90(keyGeom_);
  }
  // -----------------------------------------------------------------------------
  void Geometry::gridgen(const eckit::Configuration & config) const {
    const eckit::Configuration * conf = &config;
    Log::trace() << "Geometry::gridgen: " << keyGeom_ << std::endl;
    Log::trace() << conf << std::endl;
    soca_geo_gridgen_f90(keyGeom_, &conf);
  }
  // -----------------------------------------------------------------------------
  void Geometry::print(std::ostream & os) const {
    soca_geo_info_f90(keyGeom_);
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
