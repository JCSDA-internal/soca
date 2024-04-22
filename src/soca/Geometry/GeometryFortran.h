/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GEOMETRY_GEOMETRYFORTRAN_H_
#define SOCA_GEOMETRY_GEOMETRYFORTRAN_H_

#include "soca/Fortran.h"

#include "oops/base/Variables.h"

// Forward declarations
namespace atlas {
  namespace field {
    class FieldSetImpl;
    class FieldImpl;
  }
  namespace functionspace {
    class FunctionSpaceImpl;
  }
}
namespace eckit {
  class Configuration;
}

namespace soca {

  extern "C" {
    void soca_geo_setup_f90(F90geom &,
                            const eckit::Configuration * const &,
                            const eckit::mpi::Comm *,
                            const bool & gen);

    void soca_geo_init_atlas_f90(const F90geom &,
                      atlas::functionspace::FunctionSpaceImpl *,
                      atlas::field::FieldSetImpl *,
                      const eckit::Configuration * const &,
                      const bool & gen);
    void soca_geo_get_mesh_size_f90(const F90geom &, int &, int&);
    void soca_geo_gen_mesh_f90(
      const F90geom &,
      const int &, double[], double[], int[], int[], int[], int[],
      const int &, int[]);
    void soca_geo_delete_f90(F90geom &);
    void soca_geo_get_num_levels_f90(const F90geom &, const oops::Variables &,
                                    const size_t &, size_t[]);
    void soca_geo_iterator_dimension_f90(const F90geom &, int &);
    void soca_geo_write_f90(const F90geom &, const eckit::Configuration * const &);
  }
}  // namespace soca
#endif  // SOCA_GEOMETRY_GEOMETRYFORTRAN_H_
