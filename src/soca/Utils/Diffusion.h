/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>

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
    Diffusion( const oops::GeometryData &, const atlas::Field & scales);

    void multiply(oops::FieldSet3D &) const;
    void multiplyTL(oops::FieldSet3D &) const {}
    void multiplyAD(oops::FieldSet3D &) const {}

 private:
  atlas::Field inv_area_;
  int niterHz_;

  // stuff for preparing and storing the mesh
  const std::unique_ptr<const atlas::Mesh> createMesh(const oops::GeometryData &) const;
  const std::unique_ptr<const atlas::Mesh> mesh_;

  // Stuff for storing the calculated edge geometry
  struct EdgeGeom {
    size_t nodeA, nodeB;
    double edgeLength;
    double lengthRatio; // length of grid side / length of mesh edge
  }; 
  const std::vector<EdgeGeom> createEdgeGeom(const oops::GeometryData &) const;
  const std::vector<EdgeGeom> edgeGeom_;

  struct EdgeParam {
    std::vector<double> KhDt; 
  };
  std::vector<EdgeParam> edgeParam_;

  void multiplyHzTL(atlas::Field &) const;
  void multiplyHzAD(atlas::Field &) const;
  // void multiplyVtTL(oops::FieldSet3D &) const;
  // void multiplyVtAD(oops::FieldSet3D &) const;
  
};

// --------------------------------------------------------------------------------------

}