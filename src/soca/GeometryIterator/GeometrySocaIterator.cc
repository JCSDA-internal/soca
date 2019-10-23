/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/Configuration.h"
#include "GeometrySocaIterator.h"
#include "GeometrySocaIteratorFortran.h"
#include "Fortran.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace soca {

/*
// -----------------------------------------------------------------------------

GeometrySocaIterator::GeometrySocaIterator(const GeometrySocaIterator& iter) {
  soca_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometrySocaIterator::GeometrySocaIterator(const GeometrySoca& geom,
                                       const int & index) {
  soca_geom_iter_setup_f90(keyIter_, geom.toFortran(), index);
}


// -----------------------------------------------------------------------------

GeometrySocaIterator::~GeometrySocaIterator() {
  soca_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometrySocaIterator::operator==(const GeometrySocaIterator & other) const {
  int equals = 0;
  soca_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 1);
}

// -----------------------------------------------------------------------------

bool GeometrySocaIterator::operator!=(const GeometrySocaIterator & other) const {
  int equals = 0;
  soca_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 0);
}

// -----------------------------------------------------------------------------

eckit::geometry::Point2 GeometrySocaIterator::operator*() const {
  double lat, lon;
  soca_geom_iter_current_f90(keyIter_, lat, lon);
  return eckit::geometry::Point2(lat, lon);
}

// -----------------------------------------------------------------------------

GeometrySocaIterator& GeometrySocaIterator::operator++() {
  soca_geom_iter_next_f90(keyIter_);
  return *this;
}

// -----------------------------------------------------------------------------

void GeometrySocaIterator::print(std::ostream & os) const {
  double lat, lon;
  soca_geom_iter_current_f90(keyIter_, lat, lon);
  os << "GeometrySocaIterator, lat/lon: " << lat << " / " << lon << std::endl;
}

// -----------------------------------------------------------------------------

*/
}  // namespace soca
