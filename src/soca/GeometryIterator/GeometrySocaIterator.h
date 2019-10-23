/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef JEDI_SRC_GEOMETRYITERATOR_GEOMETRYSOCAITERATOR_H_
#define JEDI_SRC_GEOMETRYITERATOR_GEOMETRYSOCAITERATOR_H_

#include <iterator>
#include <string>

#include "eckit/geometry/Point2.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Fortran.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace soca {

class Geometry;

// -----------------------------------------------------------------------------
class GeometrySocaIterator: public std::iterator<std::forward_iterator_tag,
                                               eckit::geometry::Point2>,
                          public util::Printable,
                          private util::ObjectCounter<GeometrySocaIterator> {
 public:
  static const std::string classname() {return "soca::GeometrySocaIterator";}

  GeometrySocaIterator(const GeometrySocaIterator &);
  explicit GeometrySocaIterator(const Geometry & geom, const int & index = 1);
  ~GeometrySocaIterator();

  bool operator==(const GeometrySocaIterator &) const;
  bool operator!=(const GeometrySocaIterator &) const;
  eckit::geometry::Point2 operator*() const;
  GeometrySocaIterator& operator++();

  F90iter & toFortran() {return keyIter_;}
  const F90iter & toFortran() const {return keyIter_;}

 private:
  void print(std::ostream &) const;
  F90iter keyIter_;
};

}  // namespace soca

#endif  // JEDI_SRC_GEOMETRYITERATOR_GEOMETRYSOCAITERATOR_H_
