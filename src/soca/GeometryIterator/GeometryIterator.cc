/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"

#include "eckit/geometry/Point3.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace soca {

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator& iter)
  : geom_(iter.geom_), iIndex_(iter.iIndex_), kIndex_(iter.kIndex_) {
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry& geom,
                                   const size_t & iindex, const size_t & kindex)
 : geom_(geom), iIndex_(iindex), kIndex_(kindex) {
 }

// -----------------------------------------------------------------------------

GeometryIterator::~GeometryIterator() {}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator==(const GeometryIterator & other) const {
  return (iIndex_ == other.iIndex_ && kIndex_ == other.kIndex_);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator!=(const GeometryIterator & other) const {
  return !(*this == other );
}

// -----------------------------------------------------------------------------

eckit::geometry::Point3 GeometryIterator::operator*() const {
  ASSERT(geom_.IteratorDimension() == 2);  // modification needed for 3D
  ASSERT(iIndex_ < geom_.functionSpace().size());
  const auto & lonlat = geom_.functionSpace().lonlat();
  const auto & vLonLat = atlas::array::make_view<double,2>(lonlat);
  return eckit::geometry::Point3(vLonLat(iIndex_, 0), vLonLat(iIndex_,1), -99999);
}

// -----------------------------------------------------------------------------

double GeometryIterator::getArea() const {
  const auto & view = atlas::array::make_view<double,2>(geom_.fields().field("area"));
  return view(iIndex_, 0);
}

// -----------------------------------------------------------------------------

double GeometryIterator::getRossbyRadius() const {
  const auto & view = atlas::array::make_view<double,2>(geom_.fields().field("rossby_radius"));
  return view(iIndex_, 0);
}

// -----------------------------------------------------------------------------

GeometryIterator& GeometryIterator::operator++() {
  ASSERT(geom_.IteratorDimension() == 2);  // modification needed for 3D
  const auto & fs = geom_.functionSpace();
  if (iIndex_ >= fs.size()) {
    throw eckit::Exception("Can't go past end on geometry iterator");
  }

  const auto & ghost = atlas::array::make_view<int,1>(fs.ghost());
  do { iIndex_++; } while (iIndex_ < fs.size() && ghost(iIndex_));

  ASSERT( iIndex_ <= fs.size());
  return *this;
}

// -----------------------------------------------------------------------------

void GeometryIterator::print(std::ostream & os) const {
  const auto & p3 = **this;
  os << "GeometryIterator, lat/lon/depth: " << p3[0] << " / " << p3[1]
     << " / " << p3[2] << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace soca
