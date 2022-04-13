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
                            const eckit::mpi::Comm *);
    void soca_geo_set_atlas_lonlat_f90(const F90geom &,
                                       atlas::field::FieldSetImpl *);
    void soca_geo_set_atlas_functionspace_pointer_f90(const F90geom &,
                      atlas::functionspace::FunctionSpaceImpl *);
    void soca_geo_fill_atlas_fieldset_f90(const F90geom &,
                                          atlas::field::FieldSetImpl *);
    void soca_geo_clone_f90(F90geom &, const F90geom &);
    void soca_geo_gridgen_f90(const F90geom &);
    void soca_geo_delete_f90(F90geom &);
    void soca_geo_start_end_f90(const F90geom &, int &, int &, int &, int &,
                                int &, int &);
    void soca_geo_get_num_levels_f90(const F90geom &, const oops::Variables &,
                                    const size_t &, size_t[]);
    void soca_geo_iterator_dimension_f90(const F90geom &, int &);

    void soca_geo_gridsize_f90(const F90geom &, const char &, const bool &,
      const bool &, int &);
    void soca_geo_gridlatlon_f90(const F90geom &, const char &, const bool &,
      const bool &, const int &, double[], double[]);
    void soca_geo_getvargrid_f90(
      const F90geom &, const oops::Variables &, char &, bool &);
  }
}  // namespace soca
#endif  // SOCA_GEOMETRY_GEOMETRYFORTRAN_H_
