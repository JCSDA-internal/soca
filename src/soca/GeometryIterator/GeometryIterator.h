/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
#define SOCA_GEOMETRYITERATOR_GEOMETRYITERATOR_H_

#include <iterator>
#include <string>

#include "soca/Fortran.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  namespace geometry {
    class Point3;
  }
}
namespace soca {
  class Geometry;
}


namespace soca {
// -----------------------------------------------------------------------------
class GeometryIterator: public std::iterator<std::forward_iterator_tag,
                                               eckit::geometry::Point3>,
                          public util::Printable,
                          private util::ObjectCounter<GeometryIterator> {
 public:
  static const std::string classname() {return "soca::GeometryIterator";}

  GeometryIterator(const GeometryIterator &);
  explicit GeometryIterator(const Geometry & geom,
                            const int & iindex = 1, const int & jindex = 1,
                            const int & kindex = 1);
  ~GeometryIterator();

  bool operator==(const GeometryIterator &) const;
  bool operator!=(const GeometryIterator &) const;
  eckit::geometry::Point3 operator*() const;
  GeometryIterator& operator++();

  double getArea() const;
  double getRossbyRadius() const;

  int iteratorDimension() const;

  F90iter & toFortran() {return keyIter_;}
  const F90iter & toFortran() const {return keyIter_;}

 private:
  void print(std::ostream &) const;
  F90iter keyIter_;
};

}  // namespace soca

#endif  // SOCA_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
