/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"


#include "eckit/config/Configuration.h"
#include "eckit/config/YAMLConfiguration.h"

#include "soca/Geometry/Geometry.h"

#include "atlas/output/Gmsh.h"

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
    if (gen) {
      soca_geo_gridgen_f90(keyGeom_);
    }

    // setup the atlas functionspace
    {
      using atlas::gidx_t;
      using atlas::idx_t;

      // get the number of nodes and cells owned by this PE
      int num_nodes;
      int num_quad_elements;
      soca_geo_get_mesh_size_f90(keyGeom_, num_nodes, num_quad_elements);

      // get the mesh connectivity from the soca fortran
      std::vector<double> lons(num_nodes);
      std::vector<double> lats(num_nodes);
      std::vector<int> ghosts(num_nodes);
      std::vector<int> global_indices(num_nodes);
      std::vector<int> remote_indices(num_nodes);
      std::vector<int> partitions(num_nodes);
      const int num_quad_nodes = num_quad_elements * 4;
      std::vector<int> raw_quad_nodes(num_quad_nodes);
      soca_geo_get_mesh_f90(keyGeom_,
        num_nodes, lons.data(), lats.data(), ghosts.data(), global_indices.data(),
        remote_indices.data(), partitions.data(),
        num_quad_nodes, raw_quad_nodes.data());

      // calculate per-PE global quad numbering offset
      std::vector<int> num_elements_per_rank(comm_.size());
      comm_.allGather(num_quad_elements, num_elements_per_rank.begin(),
                      num_elements_per_rank.end());
      int global_element_index = 0;
      for (size_t i = 0; i < comm_.rank(); ++i) {
        global_element_index += num_elements_per_rank[i];
      }

      // convert some of the temporary arrays into a form atlas expects
      std::vector<gidx_t> atlas_global_indices(num_nodes);
      std::transform(global_indices.begin(), global_indices.end(), atlas_global_indices.begin(),
        [](const int index) {return atlas::gidx_t{index};});
      std::vector<idx_t> atlas_remote_indices(num_nodes);
      std::transform(remote_indices.begin(), remote_indices.end(), atlas_remote_indices.begin(),
        [](const int index) {return atlas::idx_t{index};});
      std::vector<std::array<gidx_t, 3>> tri_boundary_nodes{};  // MOM does not have triangles
      std::vector<gidx_t> tri_global_indices{};  // MOM does not have triangles
      std::vector<std::array<gidx_t, 4>> quad_boundary_nodes(num_quad_elements);
      std::vector<gidx_t> quad_global_indices(num_quad_elements);
      for (size_t quad = 0; quad < num_quad_elements; ++quad) {
        for (size_t i = 0; i < 4; ++i) {
          quad_boundary_nodes[quad][i] = raw_quad_nodes[4*quad + i];
        }
        quad_global_indices[quad] = global_element_index++;
      }

      // build the mesh!
      const atlas::idx_t remote_index_base = 1;  // 1-based indexing from Fortran
      eckit::LocalConfiguration config{};
      config.set("mpi_comm", comm_.name());
      const atlas::mesh::MeshBuilder mesh_builder{};
      atlas::Mesh mesh = mesh_builder(
        lons, lats, ghosts,
        atlas_global_indices, atlas_remote_indices, remote_index_base, partitions,
        tri_boundary_nodes, tri_global_indices,
        quad_boundary_nodes, quad_global_indices, config);
      atlas::mesh::actions::build_halo(mesh, 1);
      functionSpace_ = atlas::functionspace::NodeColumns(mesh, config);

      // save output for viewing (TODO make optional)
      atlas::output::Gmsh gmsh("out.msh",
          atlas::util::Config("coordinates", "xyz")
          | atlas::util::Config("ghost", true));  // enables viewing halos per task
      gmsh.write(mesh);
    }

    // Set ATLAS function space information in Fortran
    auto global_index = functionSpace_.global_index();
    auto ghost = functionSpace_.ghost();
    soca_geo_set_atlas_functionspace_f90(keyGeom_,
      functionSpace_.get(), global_index.get(), ghost.get(), fields_.get());
  }

  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_),
      fmsinput_(other.fmsinput_)
       {
    const int key_geo = other.keyGeom_;
    soca_geo_clone_f90(keyGeom_, key_geo);

    functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
    soca_geo_set_atlas_functionspace_f90(keyGeom_, functionSpace_.get(),
      functionSpace_.global_index().get(), functionSpace_.ghost().get(), fields_.get());
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
    // get the number of total grid points (including halo)
    int gridSizeWithHalo = functionSpace_.size();
    auto vLonlat = atlas::array::make_view<double, 2>(functionSpace_.lonlat());

    if (halo) {
      // get the latlon of all points, including halo
      lons.resize(gridSizeWithHalo);
      lats.resize(gridSizeWithHalo);
      int idx = 0;
      for (size_t i=0; i < gridSizeWithHalo; i++) {
        double lon = vLonlat(i, 0);
        double lat = vLonlat(i, 1);
        lats[idx] = lat;
        lons[idx++] = lon;
      }
    } else {
      // get the latlon of only the owned points, excluding the halo

      // count the number of owned non-ghost points (isn't there an atlas function for this??)
      auto vGhost = atlas::array::make_view<int, 1>(functionSpace_.ghost());
      int gridSize = 0;
      for (size_t i = 0; i < gridSizeWithHalo; i++) {
        if (vGhost(i) == 0) gridSize++;
      }

      // fill in the latlon
      lons.resize(gridSize);
      lats.resize(gridSize);
      int idx = 0;
      for (size_t i=0; i < gridSizeWithHalo; i++) {
        if (vGhost(i)) continue;

        double lon = vLonlat(i, 0);
        double lat = vLonlat(i, 1);
        lats[idx] = lat;
        lons[idx++] = lon;
      }
      ASSERT(idx == gridSize);
    }
  }

}  // namespace soca
