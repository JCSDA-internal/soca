/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"


#include "eckit/config/Configuration.h"
#include "eckit/config/YAMLConfiguration.h"

#include "soca/Geometry/Geometry.h"

// -----------------------------------------------------------------------------
namespace soca {

  const std::vector<char> grids{ 'h', 'u', 'v'};

  // -----------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm, const bool gen)
    : comm_(comm),
      fmsinput_(comm, conf) {

    fmsinput_.updateNameList();

    soca_geo_setup_f90(keyGeom_, &conf, &comm);
    
    // messy, fix this
    // generate the grid ONLY if being run under the gridgen application.
    // also, if true, then don't bother with the kdtree generation in the next step.
    if (gen) {
      soca_geo_gridgen_f90(keyGeom_);
      return;
    }    

    // setup the atlas functionspace
    {
      using atlas::gidx_t;
      using atlas::idx_t;
      
      int num_nodes=0;
      int num_quad_elements=0;
      soca_geo_get_mesh_size_f90(keyGeom_, num_nodes, num_quad_elements);
      num_quad_elements=0;
     
      std::vector<double> lons(num_nodes);
      std::vector<double> lats(num_nodes);
      std::vector<int> ghosts(num_nodes);
      std::vector<int> global_indices(num_nodes);
      std::vector<int> remote_indices(num_nodes);
      std::vector<int> partitions(num_nodes);

      std::vector<std::array<gidx_t, 3>> tri_boundary_nodes{};  // MOM does not have triangles
      std::vector<gidx_t> tri_global_indices{};
      std::vector<std::array<gidx_t, 4>> quad_boundary_nodes(num_quad_elements);
      std::vector<gidx_t> quad_global_indices(num_quad_elements);

      soca_geo_get_mesh_f90(keyGeom_, 
        num_nodes, lons.data(), lats.data(), ghosts.data(), global_indices.data(), 
        remote_indices.data(), partitions.data(),
        num_quad_elements);

      std::vector<gidx_t> atlas_global_indices(num_nodes);
      std::transform(global_indices.begin(), global_indices.end(), atlas_global_indices.begin(), 
        [](const int index) {return atlas::gidx_t{index};});
      std::vector<idx_t> atlas_remote_indices(num_nodes);        
      std::transform(remote_indices.begin(), remote_indices.end(), atlas_remote_indices.begin(),
        [](const int index) {return atlas::idx_t{index};});        
      
      const atlas::idx_t remote_index_base=1;  // 1-based indexing from Fortran

      eckit::LocalConfiguration config{};
      config.set("mpi_comm", comm_.name());

      std::cout << "DBG mesh_builder{}" << std::endl;
      const atlas::mesh::MeshBuilder mesh_builder{};
      std::cout << "DBG mesh_builder()" << std::endl;
      atlas::Mesh mesh = mesh_builder(
        lons, lats, ghosts, 
        atlas_global_indices, atlas_remote_indices, remote_index_base, partitions,
        tri_boundary_nodes, tri_global_indices,
        quad_boundary_nodes, quad_global_indices, config);
      std::cout << "DBG build_halo()" << std::endl;
      atlas::mesh::actions::build_halo(mesh, 1);
      std::cout << "DBG NodeColumns()" << std::endl;
      functionSpace_ = atlas::functionspace::NodeColumns(mesh, config);
      std::cout << "DBG END" << std::endl;

      ASSERT(1==2);
    }

    // // Set ATLAS lonlat and function space (with and without halos)
    // atlas::FieldSet lonlat;
    // soca_geo_lonlat_f90(keyGeom_, lonlat.get());
    // functionSpace_ = atlas::functionspace::PointCloud(lonlat->field("lonlat_inc_halos"));

    // Set ATLAS function space pointer in Fortran
    std::cout << "DBG A" << std::endl;
    soca_geo_set_atlas_functionspace_pointer_f90(keyGeom_, functionSpace_.get());

    // Fill ATLAS fieldset
    std::cout << "DBG B" << std::endl;
    // soca_geo_to_fieldset_f90(keyGeom_, fields_.get());

    std::cout << "DBG subroutine end" << std::endl;
  }
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_),
      fmsinput_(other.fmsinput_)
       {
    const int key_geo = other.keyGeom_;
    soca_geo_clone_f90(keyGeom_, key_geo);

    functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
    soca_geo_set_atlas_functionspace_pointer_f90(keyGeom_, functionSpace_.get());

    fields_ = atlas::FieldSet();
    for (int jfield = 0; jfield < other.fields_->size(); ++jfield) {
      atlas::Field atlasField = other.fields_->field(jfield);
      fields_->add(atlasField);
    }
  }
  // -----------------------------------------------------------------------------
  Geometry::~Geometry() {
    soca_geo_delete_f90(keyGeom_);
  }

  // -----------------------------------------------------------------------------
  GeometryIterator Geometry::begin() const {
    // return start of the geometry on this mpi tile
    int ist, iend, jst, jend, kst, kend, itd;
    soca_geo_start_end_f90(keyGeom_, ist, iend, jst, jend, kst, kend);
    // 3D iterator starts from 0 for surface variables
    if (IteratorDimension() == 3) kst = 0;
    return GeometryIterator(*this, ist, jst, kst);
  }
  // -----------------------------------------------------------------------------
  GeometryIterator Geometry::end() const {
    // return end of the geometry on this mpi tile
    // decided to return index out of bounds for the iterator loops to work
    return GeometryIterator(*this, -1, -1, -1);
  }
  // -----------------------------------------------------------------------------
  int Geometry::IteratorDimension() const {
    // return dimesnion of the iterator
    // if 2, iterator is over vertical columns
    // if 3, iterator is over 3D points
    int rv;
    soca_geo_iterator_dimension_f90(keyGeom_, rv);
    return rv;
  }
  // -----------------------------------------------------------------------------
  std::vector<size_t> Geometry::variableSizes(
      const oops::Variables & vars) const {
    std::vector<size_t> lvls(vars.size());
    soca_geo_get_num_levels_f90(toFortran(), vars, lvls.size(), lvls.data());
    return lvls;
  }
  // -----------------------------------------------------------------------------
  void Geometry::print(std::ostream & os) const {
    // TODO(Travis): Implement this correctly.
  }
  // -----------------------------------------------------------------------------
  void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
      const bool halo) const {
    // get the number of gridpoints
    int gridSize;
    soca_geo_gridsize_f90(keyGeom_, halo, gridSize);

    // get the lat/lon of those gridpoints
    lats.resize(gridSize);
    lons.resize(gridSize);
    soca_geo_gridlatlon_f90(keyGeom_, halo, gridSize, lats.data(), lons.data());
  }

}  // namespace soca
