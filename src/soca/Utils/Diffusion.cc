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

Diffusion::Diffusion(const oops::GeometryData & geometryData)
: geom_(geometryData),
  mesh_(createMesh(geometryData)),
  edgeGeom_(createEdgeGeom(geometryData))
{
  // calculate some other grid based constants
  //  inv_area (1/area)
  const atlas::FunctionSpace & fs = geom_.functionSpace();
  inv_area_ = fs.createField<double>();
  auto v_inv_area = atlas::array::make_view<double, 1>(inv_area_);
  const atlas::Field & area = geom_.fieldSet().field("area");
  auto v_area = atlas::array::make_view<double, 2>(area);
  for (size_t i = 0; i < fs.size(); i++) {
    v_inv_area(i) = v_area(i, 0) < 1e-6 ? 0.0 : 1.0 / v_area(i, 0);
  }
}

// --------------------------------------------------------------------------------------

void Diffusion::setScales(const atlas::FieldSet & scales) {
  const atlas::FunctionSpace & fs = geom_.functionSpace();

  scales.haloExchange();

  //-------------------------------------------------------------------------------------
  // Horizontal diffusion parameters (in units of m)
  //-------------------------------------------------------------------------------------
  // TODO(Travis) check that at least 1 of hz or vt is in fieldset
  if (scales.has("hzScales")) {
    // calculate the actual min number of iterations, and the diffusion coefficient
    const auto & hzScales = scales["hzScales"];
    const auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);

    double minItr = 0;
    khdt_.resize(edgeGeom_.size());
    khdtLevels_ = hzScales.shape(1);

    for (size_t e = 0; e < khdt_.size(); e++) {
      khdt_[e].resize(khdtLevels_);

      for (size_t level = 0; level < khdtLevels_; level++) {
        if (v_hzScales(edgeGeom_[e].nodeA, level) == 0.0 ||
            v_hzScales(edgeGeom_[e].nodeB, level) == 0.0) {
          // one of the nodes is masked out, dont diffuse with it
          khdt_[e][level] = 0.0;
        } else {
          // calculate diffusion coefficient (not taking into account the number of
          // iterations.. yet)
          const double s = (v_hzScales(edgeGeom_[e].nodeA, level) +
                            v_hzScales(edgeGeom_[e].nodeB, level)) / 2.0;
          khdt_[e][level] = s * s;

          // calculate the minimum number of iterations needed to be computationally stable
          // on this PE
          const double el = edgeGeom_[e].edgeLength;
          // TODO(Travis) check my math on this, original had a 2.0 * ... but I think 1.5 is stable
          minItr = std::max(1.5 * (s*s) / (el*el), minItr);
        }
      }
    }
    niterHz_ = round(minItr/2) * 2;  // make sure number of iterations is even
    // get the global max number of iterations
    geom_.comm().allReduceInPlace(niterHz_, eckit::mpi::Operation::MAX);

    // TODO(Travis) do some error checking to make sure niterHz_ is not too big
    std::cout << "DBG niterHz: " << niterHz_ << std::endl;

    // adjust the above calculated diffusion coefficients by the final number of
    // iterations
    for (size_t e = 0; e < khdt_.size(); e++) {
      for (size_t level = 0; level < khdtLevels_; level++) {
        khdt_[e][level] *= 1.0 / (2.0 * niterHz_);
      }
    }
  }

  //-------------------------------------------------------------------------------------
  // Vertical diffusion parameters (in units of # of levels)
  //-------------------------------------------------------------------------------------
  if (scales.has("vtScales")) {
    const auto vtScales = scales["vtScales"];
    const auto v_vtScales = atlas::array::make_view<double, 2> (vtScales);

    kvdt_ = fs.createField<double>(atlas::option::levels(vtScales.shape(1)-1));
    auto v_kvdt = atlas::array::make_view<double, 2>(kvdt_);
    double minItr = 0;
    for (size_t i = 0; i < kvdt_.shape(0); i++) {
      for (size_t level = 0; level < kvdt_.shape(1); level++) {
        if (v_vtScales(i, level) == 0.0 || v_vtScales(i, level+1) == 0.0) {
          v_kvdt(i, level) = 0.0;
        } else {
          double s = (v_vtScales(i, level) + v_vtScales(i, level+1)) / 2;
          v_kvdt(i, level) = s * s;

          // TODO(Travis) check my math on this, original had a 2.0 * ... but I think 1.5 is stable
          minItr = std::max(1.5 * s * s, minItr);
        }
      }
    }
    niterVt_ = round(minItr/2) * 2;
    geom_.comm().allReduceInPlace(niterVt_, eckit::mpi::Operation::MAX);

    // TODO(Travis) do some error checking to make sure niterVt_ is not too big
    std::cout << "DBG niterVt: " << niterVt_ << std::endl;

    // adjust the above calculated diffusion coefficients by the final number of
    // iterations
    for (size_t i = 0; i < kvdt_.shape(0); i++) {
      for (size_t lvl = 0; lvl < kvdt_.shape(1); lvl++) {
        v_kvdt(i, lvl) *= 1.0 / (2.0 * niterVt_);
      }
    }
  }
}

// --------------------------------------------------------------------------------------

const std::unique_ptr<const atlas::Mesh> Diffusion::createMesh(
    const oops::GeometryData &geometryData) const {
  // Prepare the mesh
  // TODO(Travis) talk to Francois about adding access to GeometryData::mesh_
  // Assume NodeColumns for now
  atlas::Mesh mesh;

  const auto & fs = geometryData.functionSpace();
  ASSERT(fs.type() == "NodeColumns");
  const atlas::functionspace::NodeColumns nodecolumns(fs);
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

    // NOTE there is something wrong with atlas returning a bad cell index if
    // there should only be 1 cell
    if (cellA >= mesh_->cells().size() || cellB >= mesh_->cells().size()) {
      // TODO(Travis) do this correctly
      edgeGeom.lengthRatio = 1.0;
    } else {
      const auto & centerA = atlas::Point3(centersView(cellA, 0),
                                           centersView(cellA, 1),
                                           centersView(cellA, 2));
      const auto & centerB = atlas::Point3(centersView(cellB, 0),
                                           centersView(cellB, 1),
                                           centersView(cellB, 2));
      const double sideLength = atlas::Point3::distance(centerA, centerB);
      edgeGeom.lengthRatio = edgeGeom.edgeLength < 1e-6 ? 0.0 : sideLength / edgeGeom.edgeLength;
    }
  }
  return edgeGeomVec;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiply(atlas::FieldSet &fset, Mode mode) const {
  for (atlas::Field & field : fset) {
    multiply(field);
  }
}

// --------------------------------------------------------------------------------------

void Diffusion::multiply(atlas::Field &field, Mode mode) const {
  // We do not have the true 3D diffusion operator implemented (yet?)
  // force user to use the split Hz/vt instead
  if (mode == Mode::HZVT_3D) {
    throw eckit::NotImplemented(
      "Diffusion:multiply does not have \"HZVT_3D\" mode implemented, use \"HZVT_2D_1D instead.",
      Here());
  }

  // for each field, diffuse!
  bool doVt = (mode == Mode::VT_ONLY || mode == Mode::HZVT_2D_1D)
               && field.shape(1) > 1 && niterVt_ > 0;
  bool doHz = (mode == Mode::HZ_ONLY || mode == Mode::HZVT_2D_1D) && niterHz_ > 0;

  if (doVt) multiplyVtAD(field);
  if (doHz) {
    multiplyHzAD(field);
    multiplyHzTL(field);
  }
  if (doVt) multiplyVtTL(field);
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplySqrt(atlas::FieldSet &fset, Mode mode) const {
  // We do not have the true 3D diffusion operator implemented (yet?)
  // force user to use the split Hz/vt instead
  if (mode == Mode::HZVT_3D) {
    throw eckit::NotImplemented(
      "Diffusion:multiply does not have \"HZVT_3D\" mode implemented, use \"HZVT_2D_1D instead.",
      Here());
  }

  // for each field, diffuse!
  for (atlas::Field & field : fset) {
    bool doVt = (mode == Mode::VT_ONLY || mode == Mode::HZVT_2D_1D)
                 && field.shape(1) > 1 && niterVt_ > 0;
    bool doHz = (mode == Mode::HZ_ONLY || mode == Mode::HZVT_2D_1D) && niterHz_ > 0;

    if (doHz) multiplyHzTL(field);
    if (doVt) multiplyVtTL(field);
  }
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyHzTL(atlas::Field & field) const {
  const atlas::FunctionSpace & fs = geom_.functionSpace();

  ASSERT(field.shape(0) == fs.size());
  ASSERT(khdtLevels_ == 1 || field.shape(1) <= khdtLevels_);

  const auto inv_area = atlas::array::make_view<double, 1>(inv_area_);
  auto fieldVal = atlas::array::make_view<double, 2>(field);
  std::vector<double> flux(edgeGeom_.size(), 0.0);

  for (size_t itr = 0; itr < niterHz_/2; itr++) {
    field.haloExchange();

    for (size_t level = 0; level < field.shape(1); level++) {
      // calculate diffusive flux at each edge. khdtLevels is the level to pull
      // khdt from, this code will work with either a 3D khdt field or a 2D khdt
      // field.
      const size_t khdtLevel = std::min(level, khdtLevels_-1);
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        const double dv = fieldVal(edgeGeom_[e].nodeA, level) - fieldVal(edgeGeom_[e].nodeB, level);
        flux[e] = (edgeGeom_[e].lengthRatio * dv) * khdt_[e][khdtLevel];
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
  const atlas::FunctionSpace & fs = geom_.functionSpace();

  ASSERT(field.shape(0) == fs.size());
  ASSERT(khdtLevels_ == 1 || field.shape(1) <= khdtLevels_);

  const auto & inv_area = atlas::array::make_view<double, 1>(inv_area_);
  const auto & ghost = atlas::array::make_view<int, 1>(fs.ghost());
  auto fieldVal = atlas::array::make_view<double, 2>(field);

  std::vector<double> flux(edgeGeom_.size(), 0.0);

  // init halo to zero
  for (size_t i = 0; i < fs.size(); i++) {
    if (ghost(i)) {
      for (size_t level = 0; level < field.shape(1); level++) {
        fieldVal(i, level) = 0;
      }
    }
  }

  for (size_t itr = 0; itr < niterHz_/2; itr++) {
    for (size_t level = 0; level < field.shape(1); level++) {
      // adjoint time-step the diffusion terms
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        flux[e] += inv_area(edgeGeom_[e].nodeB) * fieldVal(edgeGeom_[e].nodeB, level)
                  -inv_area(edgeGeom_[e].nodeA) * fieldVal(edgeGeom_[e].nodeA, level);
      }

      // Adjoint calculate diffusive flux at each edge. khdtLevels is the level
      // to pull khdt from, this code will work with either a 3D khdt field or a
      // 2D khdt field.
      const size_t khdtLevel = std::min(level, khdtLevels_-1);
      for (size_t e = 0; e < edgeGeom_.size(); e++) {
        fieldVal(edgeGeom_[e].nodeA, level) += edgeGeom_[e].lengthRatio
                                               * khdt_[e][khdtLevel] * flux[e];
        fieldVal(edgeGeom_[e].nodeB, level) -= edgeGeom_[e].lengthRatio
                                               * khdt_[e][khdtLevel] * flux[e];
        flux[e] = 0.0;
      }
    }

    field.adjointHaloExchange();
  }

  field.set_dirty(true);
}

// --------------------------------------------------------------------------------------

void Diffusion::multiplyVtTL(atlas::Field & field) const {
  // make sure input field is correct shape
  ASSERT(field.shape(0) == kvdt_.shape(0));
  ASSERT(field.shape(1) == kvdt_.shape(1)+1 || field.shape(1) == 1);

  field.haloExchange();
  auto v_field = atlas::array::make_view<double, 2>(field);
  const auto & v_kvdt = atlas::array::make_view<double, 2>(kvdt_);

  const int nz = field.shape(1);
  std::vector<double> flux(nz-1, 0.0);  // 1 less level than field
                                        // since this is flux between levels

  for (size_t itr = 0; itr < niterVt_ / 2; itr++) {
    for (size_t i = 0; i < field.shape(0); i++) {
      // calculate diffusive flux
      for (size_t level = 0; level < flux.size(); level++) {
        flux[level] = v_kvdt(i, level) * (v_field(i, level) - v_field(i, level+1));
      }

      // time-step diffusion terms
      for (size_t level = 0; level < nz-1; level++) {
        v_field(i, level) -= flux[level];
        v_field(i, level+1) += flux[level];
      }
    }
  }
}

// --------------------------------------------------------------------------------------
// NOTE: the vertical adjoint code should produce identical answers compared
// with the TL code above. It's explicitly coded anyway just to be safe
void Diffusion::multiplyVtAD(atlas::Field & field) const {
  // make sure input field is correct shape
  ASSERT(field.shape(0) == kvdt_.shape(0));
  ASSERT(field.shape(1) == kvdt_.shape(1)+1 || field.shape(1) == 1);

  field.haloExchange();
  auto v_field = atlas::array::make_view<double, 2>(field);
  const auto & v_kvdt = atlas::array::make_view<double, 2>(kvdt_);

  const int nz = field.shape(1);
  std::vector<double> flux(nz-1, 0.0);  // 1 less level than field
                                        // since this is flux between levels

  for (size_t itr = 0; itr < niterVt_ / 2; itr++) {
    for (size_t i = 0; i < field.shape(0); i++) {
      // adjoint time-step diffusion terms
      for (size_t level = 0; level < nz-1; level++) {
        flux[level] += v_field(i, level+1) - v_field(i, level);
      }

      // adjoint diffusive flux
      for (size_t level = 0; level < flux.size(); level++) {
        v_field(i, level) += v_kvdt(i, level) * flux[level];
        v_field(i, level+1) -= v_kvdt(i, level) * flux[level];
        flux[level] = 0.0;
      }
    }
  }
}

// --------------------------------------------------------------------------------------

}  // namespace soca
