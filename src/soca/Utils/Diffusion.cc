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
  auto area = geometryData.fieldSet().field("area");  
  auto v_area = atlas::array::make_view<double, 2>(area);
  for(size_t i = 0; i < fs.size(); i++){
    v_hfac(i) = v_area(i,0) < 1e-6 ? 0.0 : 1.0 / v_area(i,0);
  }

  // TODO clean this up
  // calculate the actual min number of iterations
  auto v_scales = atlas::array::make_view<double, 1>(scales);
  double minItr = 0;
  for(size_t e = 0; e < edgeGeom_.size(); e++){
    if (v_scales(edgeGeom_[e].nodeA) != 0.0 && v_scales(edgeGeom_[e].nodeB) != 0.0) {
      double s = (v_scales(edgeGeom_[e].nodeA) + v_scales(edgeGeom_[e].nodeB)) / 2.0;
      double el = edgeGeom_[e].edgeLength;
      // TODO check my math on this, original had a 2.0 * ... but I think this is right
      minItr = std::max(1.0 * s*s * (1.0 / (el*el)), minItr);
    }
  }
  niter_ = round(minItr/2) * 2;  // make sure number of iterations is even
  geometryData.comm().allReduceInPlace(niter_, eckit::mpi::Operation::MAX);
    
  // TODO do some error checking to make sure niter_ is not too big
  std::cout << "DBG niter: "<<niter_<<std::endl;

  // use scales to calculate KhDt
  for(size_t e = 0; e < edgeGeom_.size(); e++){
    if (v_scales(edgeGeom_[e].nodeA) != 0.0 && v_scales(edgeGeom_[e].nodeB) != 0.0) {
      double s = (v_scales(edgeGeom_[e].nodeA) + v_scales(edgeGeom_[e].nodeB)) / 2.0;
      edgeParam_[e].KhDt = s * s / (2.0 * niter_);
    } else {
      edgeParam_[e].KhDt = 0.0;
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

    // get the node indexes, make sure the lowest value index is first for
    // reasons that I might care about later, maybe
    ASSERT(edge2node.cols(i) == 2);
    const auto nodeA = edge2node(i, 0);
    const auto nodeB = edge2node(i, 1);
    ASSERT(nodeA != nodeB);

    // If both nodes are in the halo, don't bother adding this edge to the vector
    if(ghost(nodeA) && ghost(nodeB)) continue;

    // create edge in vector
    EdgeGeom &edgeGeom = edgeGeomVec.emplace_back();
    edgeGeom.nodeA = nodeA;
    edgeGeom.nodeB = nodeB;
    if(edgeGeom.nodeA > edgeGeom.nodeB) {
      std::swap(edgeGeom.nodeA, edgeGeom.nodeB);
    }

    // calculate the length of the edge
    const auto & pointA = atlas::Point3(xyz(nodeA, 0), xyz(nodeA, 1), xyz(nodeA, 2));
    const auto & pointB = atlas::Point3(xyz(nodeB, 0), xyz(nodeB, 1), xyz(nodeB, 2));
    edgeGeom.edgeLength = atlas::Point3::distance(pointA, pointB);
  
    // get cell centers, and estimate length of original model grid cell edge
    // that passes through this atlas edge.
    ASSERT(edge2cell.cols(i) == 2);
    const size_t cellA = edge2cell(i,0);
    const size_t cellB = edge2cell(i,1);
    //ASSERT(cellA >= 0); ASSERT(cellA < mesh_->cells().size());
    //ASSERT(cellB >= 0); ASSERT(cellB < mesh_->cells().size());
    // NOTE there is something wrong with atlas returning a bad cell index if there should only be 1 cell    
    if(cellA >= mesh_->cells().size() || cellB >= mesh_->cells().size()) {
      // TODO do thi correctly
      edgeGeom.lengthRatio = 1.0;
    } else {      
      const auto & centerA = atlas::Point3(centersView(cellA, 0),centersView(cellA, 1),centersView(cellA,2 ));
      const auto & centerB = atlas::Point3(centersView(cellB, 0),centersView(cellB, 1),centersView(cellB,2 ));
      double sideLength= atlas::Point3::distance(centerA,centerB);
      edgeGeom.lengthRatio = edgeGeom.edgeLength < 1e-6 ? 0.0 : sideLength / edgeGeom.edgeLength;      
    } 
  }
  return edgeGeomVec;
}

// --------------------------------------------------------------------------------------

void Diffusion::multiply(oops::FieldSet3D &fset) const {  
  const auto & node2edge = mesh_->nodes().edge_connectivity();
  auto hfac = atlas::array::make_view<double, 1>(hfac_);

  for (atlas::Field field : fset) {
    auto fieldVal = atlas::array::make_view<double, 2>(field);    
    std::vector<double> flux(edgeGeom_.size(), 0.0);

    // TODO remove this
    field.set_dirty(true);
    
    for(size_t itr=0; itr<niter_; itr++) {
      field.haloExchange();

      // calculate diffusive flux at each edge
      for(size_t e = 0; e < edgeGeom_.size(); e++){
        double dv = fieldVal(edgeGeom_[e].nodeA, 0) - fieldVal(edgeGeom_[e].nodeB, 0);
        flux[e] = (edgeGeom_[e].lengthRatio * dv) * edgeParam_[e].KhDt;
      }

      // time-step the diffusion terms
      for(size_t e = 0; e < edgeGeom_.size(); e++){
        fieldVal(edgeGeom_[e].nodeA, 0) -= hfac(edgeGeom_[e].nodeA) * flux[e];
        fieldVal(edgeGeom_[e].nodeB, 0) += hfac(edgeGeom_[e].nodeB) * flux[e];
      }
      field.set_dirty(true);
    }
  } 
}
}