/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "oops/base/GeometryData.h"

// --------------------------------------------------------------------------------------
// forward declarations
namespace atlas {
  class Field;
  class Mesh;
}
namespace oops {
  class FieldSet3D;
}

// --------------------------------------------------------------------------------------

namespace soca {

/// A pseudo-diffusion smoothing operator.
class Diffusion {
 public:
  /// The type of diffusion to perform (2D, 3D, etc..).
  enum Mode {
    HZVT_3D,  ///< A true joint 3D diffusion (not implemented, yet)
    HZVT_2D_1D,  ///< A 3D diffusion where horizontal and vertical are handled separately
    HZ_ONLY,  ///< Only horizontal diffusion is used
    VT_ONLY  ///< Only vertical diffusion is used
  };

  /// Initialize the geometry used by the diffusion operator.
  explicit Diffusion(const oops::GeometryData &);

  /// Set the horizontal and vertical scales used by the diffusion operator.
  /// @param scales contains 1 or both of the following fields:
  ///   - "hzScales" A 3D or 2D field with the horizontal length scales to use (units: meters)
  ///   - "vtScales" A 3D field with the vertical length scales to use (units: number of levels)
  void setScales(const atlas::FieldSet & scales);

  /// Perform diffusion smoothing of the input fields.
  /// The operation is guaranteed to be self-adjoint. `setScales()` must be set before
  /// this method can be called.
  void multiply(atlas::FieldSet &, Mode mode = HZVT_2D_1D) const;
  void multiply(atlas::Field &, Mode mode = HZVT_2D_1D) const;

  /// Perform the square root of the diffusion operator (e.g. the
  /// tangent-linear, with 1/2 the number of iterations compared with multiply).
  void multiplySqrt(atlas::FieldSet &, Mode mode = HZVT_2D_1D) const;

 private:
  const oops::GeometryData & geom_;

  // stuff for preparing and storing the mesh
  const std::unique_ptr<const atlas::Mesh> mesh_;
  const std::unique_ptr<const atlas::Mesh> createMesh(const oops::GeometryData &) const;

  // derived grid geometry
  atlas::Field inv_area_;
  struct EdgeGeom {
    size_t nodeA, nodeB;
    double edgeLength;
    double lengthRatio;  // length of grid side / length of mesh edge
  };
  const std::vector<EdgeGeom> edgeGeom_;
  const std::vector<EdgeGeom> createEdgeGeom(const oops::GeometryData &) const;

  // horizontal diffusion parameters
  int niterHz_ = -1;
  size_t khdtLevels_;
  std::vector<std::vector<double> > khdt_;

  // vertical diffusion parameters
  int niterVt_ = -1;
  atlas::Field kvdt_;

  // private methods where the magic happens!
  void multiplyHzTL(atlas::Field &) const;
  void multiplyHzAD(atlas::Field &) const;
  void multiplyVtTL(atlas::Field &) const;
  void multiplyVtAD(atlas::Field &) const;
};

// --------------------------------------------------------------------------------------

}  // namespace soca
