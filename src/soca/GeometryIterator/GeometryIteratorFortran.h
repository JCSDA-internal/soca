/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_
#define SOCA_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_

#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"

namespace soca {

  extern "C" {
    void soca_geom_iter_setup_f90(GeometryIterator::Ftn * &,
                                  const Geometry::Ftn * const &,
                                  const int &, const int &);
    void soca_geom_iter_clone_f90(GeometryIterator::Ftn * &,
                                  const GeometryIterator::Ftn * const &);
    void soca_geom_iter_delete_f90(GeometryIterator::Ftn * &);
    void soca_geom_iter_equals_f90(const GeometryIterator::Ftn * const &,
                                   const GeometryIterator::Ftn * const &,
                                   int &);
    void soca_geom_iter_current_f90(const GeometryIterator::Ftn * const &,
                                    double &, double &);
    void soca_geom_iter_next_f90(GeometryIterator::Ftn * &);
  }
}  // namespace soca
#endif  // SOCA_GEOMETRYITERATOR_GEOMETRYITERATORFORTRAN_H_
