/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

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
    soca_geo_setup_f90(keyGeom_, &configc);
  }
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other) {
    const int key_geo = other.keyGeom_;
    soca_geo_clone_f90(key_geo, keyGeom_);
  }
  // -----------------------------------------------------------------------------
  Geometry::~Geometry() {
    soca_geo_delete_f90(keyGeom_);
  }

  // -----------------------------------------------------------------------------          /**/   
  void Geometry::print(std::ostream & os) const {
    soca_geo_info_f90(keyGeom_);

  }
  // -----------------------------------------------------------------------------
}  // namespace soca
