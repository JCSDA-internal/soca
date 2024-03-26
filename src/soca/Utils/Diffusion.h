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

class Diffusion {
 public:
    Diffusion(const oops::GeometryData &,
              const atlas::Field & hzScales,
              const atlas::Field & vtScales);

    void multiply(oops::FieldSet3D &) const;
    // void multiplyTL(oops::FieldSet3D &) const {}
    // void multiplyAD(oops::FieldSet3D &) const {}

 private:
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
  int niterHz_;
  std::vector<std::vector<double> > khdt_;

  // vertical diffusion parameters
  int niterVt_;
  atlas::Field kvdt_;

  // private methods where the magic happens!
  void multiplyHzTL(atlas::Field &) const;
  void multiplyHzAD(atlas::Field &) const;
  void multiplyVtTL(atlas::Field &) const;
  // void multiplyVtAD(oops::FieldSet3D &) const;
};

// --------------------------------------------------------------------------------------

}  // namespace soca
