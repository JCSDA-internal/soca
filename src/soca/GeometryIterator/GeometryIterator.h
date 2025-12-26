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
/// \brief Iterator over the SOCA grid points
/// \details
/// This iterator allows traversal of the 3D grid points in the SOCA Geometry object.
/// It iterates over non-ghost horizontal points (indexed by i) and, if 3D iterator is
/// requested in Geometry, over vertical levels (indexed by k).
///
/// The iterator provides methods to:
/// - Access 3D coordinates (longitude, latitude, vertical coordinate) via operator*()
/// - Retrieve 2D field values at the current position via getFieldValue()
/// - Get current i and k indices directly
///
/// Vertical coordinate values are retrieved from the "vert_coord" field in the
/// Geometry object.
class GeometryIterator:   public util::Printable,
                          private util::ObjectCounter<GeometryIterator> {
 public:
  /// Standard iterator type definitions to comply with STL requirements
  typedef std::forward_iterator_tag iterator_category;
  typedef eckit::geometry::Point3 value_type;
  typedef ptrdiff_t difference_type;
  typedef eckit::geometry::Point3 & reference;
  typedef eckit::geometry::Point3 * pointer;

  /// \brief Returns the class name for logging and identification
  static const std::string classname() {return "soca::GeometryIterator";}

  /// \brief Copy constructor
  /// \param iter Iterator to copy from
  GeometryIterator(const GeometryIterator & iter);

  /// \brief Constructor with explicit position
  /// \param geom Reference to the SOCA Geometry object
  /// \param iindex Horizontal grid point index (0-based)
  /// \param kindex Vertical level index (0-based, defaults to 0 for 2D iterators)
  explicit GeometryIterator(const Geometry & geom,
                            const size_t & iindex,
                            const size_t & kindex = 0);

  /// \brief Destructor
  ~GeometryIterator() = default;

  /// \brief Equality comparison operator
  /// \param other Iterator to compare with
  /// \return True if both iterators point to the same position
  bool operator==(const GeometryIterator & other) const;

  /// \brief Inequality comparison operator
  /// \param other Iterator to compare with
  /// \return True if iterators point to different positions
  bool operator!=(const GeometryIterator & other) const;

  /// \brief Dereference operator, returns 3D point coordinates
  /// \return Point3 object containing (lon, lat, vert_coord) coordinates
  eckit::geometry::Point3 operator*() const;

  /// \brief Prefix increment operator, advances to next grid point
  /// \return Reference to this iterator after increment
  GeometryIterator& operator++();

  /// \brief Retrieves 2D field value at current position
  /// \param fieldname Name of the field to retrieve
  /// \return Value of the specified field at current position
  double getFieldValue(const std::string & fieldname) const;

  /// \brief Get current horizontal index
  /// \return Current i-index (horizontal position)
  size_t i() const {return iIndex_;}

  /// \brief Get current vertical level index
  /// \return Current k-index (vertical level)
  size_t k() const {return kIndex_;}

 private:
  /// \brief Implementation of Printable interface
  /// \param os Output stream
  void print(std::ostream & os) const;

  const Geometry & geom_;  ///< Reference to SOCA Geometry object
  size_t iIndex_;          ///< Current horizontal grid point index
  size_t kIndex_;          ///< Current vertical level index
  size_t klev_;            ///< Total number of vertical levels
};

}  // namespace soca

