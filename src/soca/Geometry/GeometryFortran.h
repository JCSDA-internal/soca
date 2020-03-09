/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GEOMETRY_GEOMETRYFORTRAN_H_
#define SOCA_GEOMETRY_GEOMETRYFORTRAN_H_

#include "soca/Geometry/Geometry.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_geo_setup_f90(Geometry::Ftn * &,
                            const eckit::Configuration * const *);
    void soca_geo_clone_f90(Geometry::Ftn * &, const Geometry::Ftn * const &);
    void soca_geo_gridgen_f90(const Geometry::Ftn * const &);
    void soca_geo_delete_f90(Geometry::Ftn * &);
    void soca_geo_start_end_f90(const Geometry::Ftn * const &,
                                int &, int &, int &, int &);
  }
}  // namespace soca
#endif  // SOCA_GEOMETRY_GEOMETRYFORTRAN_H_
