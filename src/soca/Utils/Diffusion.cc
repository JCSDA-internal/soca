#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildXYZField.h"

#include "oops/base/FieldSet3D.h"

#include "soca/Utils/Diffusion.h"

namespace soca {

// --------------------------------------------------------------------------------------

Diffusion::Diffusion(const oops::GeometryData & geometryData, const atlas::Field & scales)
: mesh_(createMesh(geometryData)),
  edgeGeom_(createEdgeGeom(geometryData)),
  edgeParam_(edgeGeom_.size())
{
  
  // calc hfac (1/area)
  const atlas::FunctionSpace & fs = geometryData.functionSpace();
  hfac_ = fs.createField<double>();
  auto v_hfac = atlas::array::make_view<double, 1>(hfac_);
  auto v_area = atlas::array::make_view<double, 2>(geometryData.fieldSet().field("area"));
  for(size_t i = 0; i < fs.size(); i++){
    v_hfac(i) = v_area(i,0) > 1e-9 ? 1.0 / v_area(i,0) : 0.0;
  }

  // TODO calculate the actual min number of iterations
  niter_ = 200;

  // use scales to calculate KhDt
  auto v_scales = atlas::array::make_view<double, 1>(scales);
  for(size_t e = 0; e < edgeGeom_.size(); e++){
    if (v_scales(edgeGeom_[e].nodeA) != 0.0 && v_scales(edgeGeom_[e].nodeB) != 0.0) {
      double s = (v_scales(edgeGeom_[e].nodeA) + v_scales(edgeGeom_[e].nodeB)) / 2.0;
      edgeParam_[e].KhDt = s * s / (2.0 * niter_);
    }
  }
  // // TODO grab some fields
  // // area, gmask?
  // for (auto & f : geometryData.fieldSet()) {
  //   std::cout << "DBG " << f.name() << std::endl;
  // }

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

const std::vector<Diffusion::EdgeGeom> Diffusion::createEdgeGeom(const oops::GeometryData &geometryData) const {  
  // process the geometry. In short, for each edge of the mesh, we need to
  // calculate the length of the edge (easy) and the length of the grid box side
  // that crosses this edge (not as easy).
  std::vector<EdgeGeom> edgeGeomVec;
 
  // get the fields we'll need later
  const auto ghost = atlas::array::make_view<int,1>(mesh_->nodes().ghost());
  const auto xyz = atlas::array::make_view<double,2>(mesh_->nodes().field("xyz"));
  const auto centersView = atlas::array::make_view<double, 2>(mesh_->cells().field("centre"));

  // calculate grid metrics on the edges
  const auto & edge2node = mesh_->edges().node_connectivity();
  const auto & edge2cell = mesh_->edges().cell_connectivity();

  edgeGeomVec.reserve(mesh_->edges().size());
  for (size_t i = 0; i < mesh_->edges().size(); i++) {
    EdgeGeom &edgeGeom = edgeGeomVec.emplace_back();

    // get the node indexes, make sure the lowest value index is first for
    // reasons that I might care about later, maybe
    ASSERT(edge2node.cols(i) == 2);
    const auto & nodeA = edge2node(i, 0); edgeGeom.nodeA = nodeA;
    const auto & nodeB = edge2node(i, 1); edgeGeom.nodeB = nodeB;
    ASSERT(edgeGeom.nodeA != edgeGeom.nodeB);
    if(edgeGeom.nodeA > edgeGeom.nodeB) {
      std::swap(edgeGeom.nodeA, edgeGeom.nodeB);
    }

    // make sure both nodes aren't in the halo
    edgeGeom.valid = !(ghost(edgeGeom.nodeA) && ghost(edgeGeom.nodeB));    
    if (!edgeGeom.valid) continue;

    // calculate the length of the edge
    const auto & pointA = atlas::Point3(xyz(nodeA, 0), xyz(nodeA, 1), xyz(nodeA, 2));
    const auto & pointB = atlas::Point3(xyz(nodeB, 0), xyz(nodeB, 1), xyz(nodeB, 2));
    edgeGeom.edgeLength = atlas::Point3::distance(pointA, pointB);
  
    // get cell centers, and estimate length of original model grid cell edge
    // that passes through this atlas edge.
    ASSERT(edge2cell.cols(i) == 2);
    const size_t cellA = edge2cell(i,0);
    const size_t cellB = edge2cell(i,1);
    ASSERT(cellA >= 0);
    ASSERT(cellB >= 0);
    const auto & centerA = atlas::Point3(centersView(cellA, 0),centersView(cellA, 1),centersView(cellA,2 ));
    const auto & centerB = atlas::Point3(centersView(cellB, 0),centersView(cellB, 1),centersView(cellB,2 ));
    edgeGeom.sideLength= atlas::Point3::distance(centerA,centerB);

    edgeGeom.dsde = edgeGeom.edgeLength < 1e-6 ? 0.0 : edgeGeom.sideLength / edgeGeom.edgeLength;
  }
  return edgeGeomVec;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiply(oops::FieldSet3D &fset) const {

  // get edge connectivity information
  const auto & node2edge = mesh_->nodes().edge_connectivity();

  auto hfac = atlas::array::make_view<double, 1>(hfac_);

  for (atlas::Field field : fset) {

    auto fieldVal = atlas::array::make_view<double, 2>(field);

    // TEMP
    for(size_t i =0; i < field.shape(0); i++){
      fieldVal(i,1) =   fieldVal(i,0);
    }
    
    // calculate diffusive flux at each edge
    std::vector<double> flux(edgeGeom_.size(), 0.0);

    field.haloExchange();
    for(size_t itr=0; itr<niter_; itr++) {
      for(size_t e = 0; e < edgeGeom_.size(); e++){
        double dv = fieldVal(edgeGeom_[e].nodeA, 0) - fieldVal(edgeGeom_[e].nodeB, 0);
        flux[e] = edgeGeom_[e].dsde * edgeParam_[e].KhDt * dv;
      }

      // // TEMP
      // for(size_t i =0; i < field.shape(0); i++){
      //   fieldVal(i,0) = 0.0;
      // }

      for(size_t e = 0; e < edgeGeom_.size(); e++){
        if (! edgeGeom_[e].valid) continue;

        ASSERT(edgeGeom_[e].nodeA >= 0);
        ASSERT(edgeGeom_[e].nodeB >= 0);
        ASSERT(edgeGeom_[e].nodeA < field.shape(0));
        ASSERT(edgeGeom_[e].nodeB < field.shape(0));

        fieldVal(edgeGeom_[e].nodeA, 0) -= hfac(edgeGeom_[e].nodeA) * flux[e];
        fieldVal(edgeGeom_[e].nodeB, 0) += hfac(edgeGeom_[e].nodeB) * flux[e];
        //fieldVal(edgeGeom_[e].nodeA, 0) = edgeParam_[e].KhDt;
      }

      field.haloExchange();
    }
  }
  

}

}