/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/Configuration.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/GeometryIterator/GeometryIteratorFortran.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace soca {


// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator& iter) {
  soca_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry& geom,
                                       const int & iindex, const int & jindex) {
  soca_geom_iter_setup_f90(keyIter_, geom.toFortran(), iindex, jindex);
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

eckit::geometry::Point2 GeometryIterator::operator*() const {
  double lat, lon;
  soca_geom_iter_current_f90(keyIter_, lat, lon);
  return eckit::geometry::Point2(lon, lat);
}

// -----------------------------------------------------------------------------

GeometryIterator& GeometryIterator::operator++() {
  soca_geom_iter_next_f90(keyIter_);
  return *this;
}

// -----------------------------------------------------------------------------

void GeometryIterator::print(std::ostream & os) const {
  double lat, lon;
  soca_geom_iter_current_f90(keyIter_, lat, lon);
  os << "GeometryIterator, lat/lon: " << lat << " / " << lon << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace soca
