/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "atlas/functionspace.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/output/Gmsh.h"

#include "eckit/config/Configuration.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Utils/readNcAndInterp.h"

// -----------------------------------------------------------------------------
namespace soca {

  // -----------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm, const bool gen)
    : comm_(comm),
      fmsinput_(comm, conf) {

    fmsinput_.updateNameList();

    // Create the grid decomposition from MOM6.
    // Also either use MOM6 to generate the grid fields, or read precomputed values in.
    soca_geo_setup_f90(keyGeom_, &conf, &comm, gen);

    // setup the atlas functionSpace
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
      soca_geo_gen_mesh_f90(keyGeom_,
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

      // optionally save output for viewing with gmsh
      if (conf.getBool("gmsh save", false)) {
        std::string filename = conf.getString("gmsh filename", "out.msh");
        atlas::output::Gmsh gmsh(filename,
            atlas::util::Config("coordinates", "xyz")
            | atlas::util::Config("ghost", true));  // enables viewing halos per task
        gmsh.write(mesh);
      }
    }

    // Set ATLAS function space in Fortran, and either generate some of the atlas fields
    // if "gen" is set, otherwise those fields will be read in from the gridspec file
    soca_geo_init_atlas_f90(keyGeom_, functionSpace_.get(), fields_.get(), &conf, gen);

    // fill in parts of the fieldset from the C++ side here
    if (gen) {
      // calculate rossby radius
      std::string rossbyFile = conf.getString("rossby file");
      auto results = readNcAndInterp(rossbyFile, {"rossby_radius"}, functionSpace_);
      fields_.add(results.field(0));
    }

    // write output
    if (gen) {
      soca_geo_write_f90(keyGeom_, &conf);
    }
  }

  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_), fmsinput_(other.fmsinput_) {
    throw eckit::Exception("Geometry copy constructor is not implemented");
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

}  // namespace soca
