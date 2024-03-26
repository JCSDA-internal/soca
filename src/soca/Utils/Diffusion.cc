/*
 * (C) Copyright 2024-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <algorithm>
#include <utility>

#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildXYZField.h"

#include "oops/base/FieldSet3D.h"

#include "soca/Utils/Diffusion.h"

namespace soca {

// --------------------------------------------------------------------------------------

Diffusion::Diffusion(const oops::GeometryData & geometryData,
                     const atlas::Field & hzScales,
                     const atlas::Field & vtScales)
: mesh_(createMesh(geometryData)),
  edgeGeom_(createEdgeGeom(geometryData))
{
  // calculate some other grid based constants
  //  inv_area (1/area)
  const atlas::FunctionSpace & fs = geometryData.functionSpace();
  inv_area_ = fs.createField<double>();
  auto v_inv_area = atlas::array::make_view<double, 1>(inv_area_);
  auto area = geometryData.fieldSet().field("area");
  auto v_area = atlas::array::make_view<double, 2>(area);
  for (size_t i = 0; i < fs.size(); i++) {
    v_inv_area(i) = v_area(i, 0) < 1e-6 ? 0.0 : 1.0 / v_area(i, 0);
  }

  //-------------------------------------------------------------------------------------
  // Horizontal diffusion parameters (in units of m)
  //-------------------------------------------------------------------------------------

  // calculate the actual min number of iterations, and the diffusion coefficient
  auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);
  double minItr = 0;
  khdt_.resize(edgeGeom_.size());
  for (size_t e = 0; e < khdt_.size(); e++) {
    khdt_[e].resize(hzScales.shape(1));

    for (size_t level = 0; level < hzScales.shape(1); level++) {
      if (v_hzScales(edgeGeom_[e].nodeA, level) == 0.0 ||
          v_hzScales(edgeGeom_[e].nodeB, level) == 0.0) {
        // one of the nodes must be masked out, dont diffuse with it
        khdt_[e][level] = 0.0;
      } else {
        // calculate diffusion coefficient (not taking into account the number of
        // iterations.. yet)
        double s = (v_hzScales(edgeGeom_[e].nodeA, level) +
                    v_hzScales(edgeGeom_[e].nodeB, level)) / 2.0;
        khdt_[e][level] = s * s;

        // calculate the minimum number of iterations needed to be computationally stable
        // on this PE
        double el = edgeGeom_[e].edgeLength;
        // TODO check my math on this, original had a 2.0 * ... but I think this is right
        minItr = std::max(1.0 * (s*s) / (el*el), minItr);
      }
    }
  }
  niterHz_ = round(minItr/2) * 2;  // make sure number of iterations is even
  // get the global max number of iterations
  geometryData.comm().allReduceInPlace(niterHz_, eckit::mpi::Operation::MAX);

  // TODO do some error checking to make sure niterHz_ is not too big
  std::cout << "DBG niterHz: " << niterHz_ << std::endl;

  // adjust the above calculated diffusion coefficients by the final number of
  // iterations
  for (size_t e = 0; e < khdt_.size(); e++) {
    for (size_t level = 0; level < khdt_[e].size(); level++) {
      khdt_[e][level] *= 1.0 / (2.0 * niterHz_);
    }
  }

  //-------------------------------------------------------------------------------------
  // Vertical diffusion parameters (in units of # of levels)
  //-------------------------------------------------------------------------------------
  kvdt_ = fs.createField<double>(atlas::option::levels(vtScales.shape(1)-1));
  auto v_kvdt = atlas::array::make_view<double, 2>(kvdt_);
  auto v_vtScales = atlas::array::make_view<double, 2> (vtScales);
  minItr = 0;
  for (size_t i = 0; i < kvdt_.shape(0); i++) {
    for (size_t level = 0; level < kvdt_.shape(1); level++) {
      if (v_vtScales(i, level) == 0.0 || v_vtScales(i, level+1) == 0.0) {
        v_kvdt(i, level) = 0.0;
      } else {
        double s = (v_vtScales(i, level) + v_vtScales(i, level+1)) / 2;
        v_kvdt(i, level) = s * s;

        minItr = std::max(1.0 * s * s, minItr);  // check math, should this be 2.0?
      }
    }
  }
  niterVt_ = round(minItr/2) * 2;
  geometryData.comm().allReduceInPlace(niterVt_, eckit::mpi::Operation::MAX);

  std::cout << "DBG niterVt: " << niterVt_ << std::endl;

  // adjust the above calculated diffusion coefficients by the final number of
  // iterations
  for (size_t i = 0; i < kvdt_.shape(0); i++) {
    for (size_t lvl = 0; lvl < kvdt_.shape(1); lvl++) {
      v_kvdt(i, lvl) *= 1.0 / (2.0 * niterVt_);
    }
  }
}

// --------------------------------------------------------------------------------------

const std::unique_ptr<const atlas::Mesh> Diffusion::createMesh(const oops::GeometryData &geometryData) const {
  // Prepare the mesh
  // TODO talk to Francois about adding access to GeometryData::mesh_
  // Assume NodeColumns for now
  atlas::Mesh mesh;

  auto fs = geometryData.functionSpace();
  ASSERT(fs.type() == "NodeColumns");
  atlas::functionspace::NodeColumns nodecolumns(fs);
  mesh = nodecolumns.mesh();
  atlas::mesh::actions::build_edges(mesh);
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  atlas::mesh::actions::BuildXYZField()(mesh);
  atlas::mesh::actions::BuildCellCentres()(mesh);
  ASSERT(mesh.nodes().has_field("xyz"));
  ASSERT(mesh.cells().has_field("centre"));

  return std::make_unique<const atlas::Mesh>(mesh);
}

// --------------------------------------------------------------------------------------

const std::vector<Diffusion::EdgeGeom> Diffusion::createEdgeGeom(
    const oops::GeometryData &geometryData) const {
  // process the geometry. In short, for each edge of the mesh, we need to
  // calculate the length of the edge (easy) and the length of the grid box side
  // that crosses this edge (not as easy).
  std::vector<EdgeGeom> edgeGeomVec;

  // get the fields we'll need later
  const auto ghost = atlas::array::make_view<int, 1>(mesh_->nodes().ghost());
  const auto xyz = atlas::array::make_view<double, 2>(mesh_->nodes().field("xyz"));
  const auto centersView = atlas::array::make_view<double, 2>(mesh_->cells().field("centre"));

  // calculate grid metrics on the edges
  const auto & edge2node = mesh_->edges().node_connectivity();
  const auto & edge2cell = mesh_->edges().cell_connectivity();

  edgeGeomVec.reserve(mesh_->edges().size());
  for (size_t i = 0; i < mesh_->edges().size(); i++) {
    // get the node indexes, make sure the lowest value index is first for
    // reasons that I might care about later, maybe
    ASSERT(edge2node.cols(i) == 2);
    const auto nodeA = edge2node(i, 0);
    const auto nodeB = edge2node(i, 1);
    ASSERT(nodeA != nodeB);

    // If both nodes are in the halo, don't bother adding this edge to the vector
    if (ghost(nodeA) && ghost(nodeB)) continue;

    // create edge in vector
    EdgeGeom &edgeGeom = edgeGeomVec.emplace_back();
    edgeGeom.nodeA = nodeA;
    edgeGeom.nodeB = nodeB;
    if (edgeGeom.nodeA > edgeGeom.nodeB) {
      std::swap(edgeGeom.nodeA, edgeGeom.nodeB);
    }

    // calculate the length of the edge
    const auto & pointA = atlas::Point3(xyz(nodeA, 0), xyz(nodeA, 1), xyz(nodeA, 2));
    const auto & pointB = atlas::Point3(xyz(nodeB, 0), xyz(nodeB, 1), xyz(nodeB, 2));
    edgeGeom.edgeLength = atlas::Point3::distance(pointA, pointB);

    // get cell centers, and estimate length of original model grid cell edge
    // that passes through this atlas edge.
    ASSERT(edge2cell.cols(i) == 2);
    const size_t cellA = edge2cell(i, 0);
    const size_t cellB = edge2cell(i, 1);
    // ASSERT(cellA >= 0); ASSERT(cellA < mesh_->cells().size());
    // ASSERT(cellB >= 0); ASSERT(cellB < mesh_->cells().size());
    // NOTE there is something wrong with atlas returning a bad cell index if there should only be 1 cell
    if (cellA >= mesh_->cells().size() || cellB >= mesh_->cells().size()) {
      // TODO do this correctly
      edgeGeom.lengthRatio = 1.0;
    } else {
      const auto & centerA = atlas::Point3(centersView(cellA, 0), centersView(cellA, 1), centersView(cellA, 2));
      const auto & centerB = atlas::Point3(centersView(cellB, 0), centersView(cellB, 1), centersView(cellB, 2));
      double sideLength = atlas::Point3::distance(centerA, centerB);
      edgeGeom.lengthRatio = edgeGeom.edgeLength < 1e-6 ? 0.0 : sideLength / edgeGeom.edgeLength;
    }
  }
  return edgeGeomVec;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiply(oops::FieldSet3D &fset) const {
  for (atlas::Field field : fset) {
    // multiplyHzAD(field);
    multiplyVtTL(field);
    multiplyHzTL(field);
    multiplyHzTL(field);
    multiplyVtTL(field);
  }
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyHzTL(atlas::Field & field) const {
  auto inv_area = atlas::array::make_view<double, 1>(inv_area_);
  auto fieldVal = atlas::array::make_view<double, 2>(field);
  std::vector<double> flux(edgeGeom_.size(), 0.0);

  // TODO remove this
  field.set_dirty(true);

  for (size_t itr = 0; itr < niterHz_/2; itr++) {
    field.haloExchange();

    for (size_t level = 0; level < field.shape(1); level++) {
      // TODO, also handle case where 2D scales are applied to 3D field??
      // TODO, invert the order of loops? go over each edge first, then the
      // level within each edge?

      // calculate diffusive flux at each edge
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        double dv = fieldVal(edgeGeom_[e].nodeA, level) - fieldVal(edgeGeom_[e].nodeB, level);
        flux[e] = (edgeGeom_[e].lengthRatio * dv) * khdt_[e][level];
      }

      // time-step the diffusion terms
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        fieldVal(edgeGeom_[e].nodeA, level) -= inv_area(edgeGeom_[e].nodeA) * flux[e];
        fieldVal(edgeGeom_[e].nodeB, level) += inv_area(edgeGeom_[e].nodeB) * flux[e];
      }
    }
    field.set_dirty(true);
  }
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyHzAD(atlas::Field & field) const {
  ASSERT(1 == 2);
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyVtTL(atlas::Field & field) const {
  auto v_field = atlas::array::make_view<double, 2>(field);
  auto v_kvdt = atlas::array::make_view<double, 2>(kvdt_);

  const int nz = field.shape(1);
  std::vector<double> flux(nz-1, 0.0);  // 1 less level than field, since this is flux between levels

  field.set_dirty(true);
  for (size_t itr = 0; itr < niterVt_ / 2; itr++) {
    field.haloExchange();

    for (size_t i = 0; i < field.shape(0); i++) {
      // calculate diffusive flux
      for (size_t level = 0; level < flux.size(); level++) {
        flux[level] = v_kvdt(i, level) * (v_field(i, level) - v_field(i, level+1)) / 2;
      }

      // time-step diffusion terms
      for (size_t level = 0; level < nz-1; level++) {
        v_field(i, level) -= flux[level];
        v_field(i, level+1) += flux[level];
      }
    }
    field.set_dirty(true);
  }
}

// --------------------------------------------------------------------------------------

}  // namespace soca
