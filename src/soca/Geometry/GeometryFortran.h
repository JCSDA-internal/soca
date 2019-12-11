/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GEOMETRY_GEOMETRYFORTRAN_H_
#define SOCA_GEOMETRY_GEOMETRYFORTRAN_H_

#include "soca/Fortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
    void soca_geo_clone_f90(const F90geom &, F90geom &);
    void soca_geo_gridgen_f90(const F90geom &);
    void soca_geo_delete_f90(F90geom &);
    void soca_geo_start_end_f90(const F90geom &, int &, int &, int &, int &);
    void soca_geo_global_grid_size_f90(const F90geom &, int &, int &, int &);
  }
}  // namespace soca
#endif  // SOCA_GEOMETRY_GEOMETRYFORTRAN_H_
