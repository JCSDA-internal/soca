/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/GeometryIterator/GeometryIteratorFortran.h"

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace soca {


// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator& iter) {
  soca_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry& geom,
                                   const int & iindex, const int & jindex,
                                   const int & kindex) {
  soca_geom_iter_setup_f90(keyIter_, geom.toFortran(), iindex, jindex, kindex);
}


// -----------------------------------------------------------------------------

GeometryIterator::~GeometryIterator() {
  soca_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator==(const GeometryIterator & other) const {
  int equals = 0;
  soca_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 1);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator!=(const GeometryIterator & other) const {
  int equals = 0;
  soca_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 0);
}

// -----------------------------------------------------------------------------

eckit::geometry::Point3 GeometryIterator::operator*() const {
  double lat, lon, dep;
  soca_geom_iter_current_f90(keyIter_, lon, lat, dep);
  return eckit::geometry::Point3(lon, lat, dep);
}

// -----------------------------------------------------------------------------

double GeometryIterator::getArea() const {
  double val;
  soca_geom_iter_get_area_f90(keyIter_, val);
  return val;
}

// -----------------------------------------------------------------------------

double GeometryIterator::getRossbyRadius() const {
  double val;
  soca_geom_iter_get_rossby_f90(keyIter_, val);
  return val;
}
// -----------------------------------------------------------------------------

GeometryIterator& GeometryIterator::operator++() {
  soca_geom_iter_next_f90(keyIter_);
  return *this;
}
// -----------------------------------------------------------------------------

int GeometryIterator::iteratorDimension() const {
  int dimension;
  soca_geom_iter_dimension_f90(keyIter_, dimension);
  return dimension;
}

// -----------------------------------------------------------------------------

void GeometryIterator::print(std::ostream & os) const {
  double lat, lon, dep;
  soca_geom_iter_current_f90(keyIter_, lon, lat, dep);
  os << "GeometryIterator, lat/lon/depth: " << lat << " / " << lon
     << " / " << dep << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace soca
