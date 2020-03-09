/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "soca/Geometry/Geometry.h"
#include "soca/Geometry/GeometryFortran.h"
#include "soca/GeometryIterator/GeometryIterator.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm), atmconf_(conf), initatm_(initAtm(conf)) {
    const eckit::Configuration * configc = &conf;
    soca_geo_setup_f90(ftn_, &configc);
  }
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_),
      atmconf_(other.atmconf_),
      initatm_(initAtm(other.atmconf_)) {
    soca_geo_clone_f90(ftn_, other.ftn_);
  }
  // -----------------------------------------------------------------------------
  Geometry::~Geometry() {
    soca_geo_delete_f90(ftn_);
  }
  // -----------------------------------------------------------------------------
  void Geometry::gridgen(const eckit::Configuration & config) const {
    Log::trace() << "Geometry::gridgen: " << std::endl;
    Log::trace() << config << std::endl;
    soca_geo_gridgen_f90(ftn_);
  }
  // -----------------------------------------------------------------------------
  GeometryIterator Geometry::begin() const {
    // return start of the geometry on this mpi tile
    int ist, iend, jst, jend;
    soca_geo_start_end_f90(ftn_, ist, iend, jst, jend);
    return GeometryIterator(*this, ist, jst);
  }
  // -----------------------------------------------------------------------------
  GeometryIterator Geometry::end() const {
    // return end of the geometry on this mpi tile
    // decided to return index out of bounds for the iterator loops to work
    return GeometryIterator(*this, -1, -1);
  }
  // -----------------------------------------------------------------------------
  void Geometry::print(std::ostream & os) const {
    // TODO(Travis): Implement this correctly.
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
