/*
 * (C) Copyright 2019-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iterator>
#include <string>

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
class GeometryIterator:   public util::Printable,
                          private util::ObjectCounter<GeometryIterator> {
 public:
  typedef std::forward_iterator_tag iterator_category;
  typedef eckit::geometry::Point3 value_type;
  typedef ptrdiff_t difference_type;
  typedef eckit::geometry::Point3 & reference;
  typedef eckit::geometry::Point3 * pointer;

  static const std::string classname() {return "soca::GeometryIterator";}

  GeometryIterator(const GeometryIterator &);
  explicit GeometryIterator(const Geometry & geom,
                            const size_t & iindex,
                            const size_t & kindex = -1);
  ~GeometryIterator();

  bool operator==(const GeometryIterator &) const;
  bool operator!=(const GeometryIterator &) const;
  eckit::geometry::Point3 operator*() const;
  GeometryIterator& operator++();

  // TODO(Travis) generalize this to get any geom field
  double getArea() const;
  double getRossbyRadius() const;

  const size_t i() const {return iIndex_;}
  const size_t k() const {return kIndex_;}

 private:
  void print(std::ostream &) const;
  const Geometry & geom_;
  size_t iIndex_;
  size_t kIndex_;
};

}  // namespace soca

