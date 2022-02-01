/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "eckit/config/YAMLConfiguration.h"

#include "soca/Geometry/Geometry.h"

// -----------------------------------------------------------------------------
namespace soca {
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm),
      fmsinput_(comm, conf) {

    fmsinput_.updateNameList();

    soca_geo_setup_f90(keyGeom_, &conf, &comm);

    // Set ATLAS lon/lat field
    atlas::FieldSet atlasFieldSet;
    soca_geo_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet.get());
    atlas::Field atlasField = atlasFieldSet.field("lonlat");

    // Create ATLAS function space
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));
//    atlasFunctionSpaceHalo_.reset(new atlas::functionspace::PointCloud(atlasField,
//      atlas::option::halo(1))); TBD: no halo available for the PointCloud function space
    atlasFunctionSpaceHalo_.reset(new atlas::functionspace::PointCloud(atlasField));

    // Set ATLAS function space pointer in Fortran
    soca_geo_set_atlas_functionspace_pointer_f90(keyGeom_,
      atlasFunctionSpace_->get());

    // Fill ATLAS fieldset
    atlasFieldSet_.reset(new atlas::FieldSet());
    soca_geo_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());
  }
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_),
      fmsinput_(other.fmsinput_) {
    const int key_geo = other.keyGeom_;
    soca_geo_clone_f90(keyGeom_, key_geo);
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                              other.atlasFunctionSpace_->lonlat()));
//    atlasFunctionSpaceHalo_.reset(new atlas::functionspace::PointCloud(
//                              other.atlasFunctionSpaceHalo_->lonlat(),atlas::option::halo(1)));
    atlasFunctionSpaceHalo_.reset(new atlas::functionspace::PointCloud(
                            other.atlasFunctionSpaceHalo_->lonlat()));
    soca_geo_set_atlas_functionspace_pointer_f90(keyGeom_,
      atlasFunctionSpace_->get());
    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
      atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
      atlasFieldSet_->add(atlasField);
    }
  }
  // -----------------------------------------------------------------------------
  Geometry::~Geometry() {
    soca_geo_delete_f90(keyGeom_);
  }
  // -----------------------------------------------------------------------------
  void Geometry::gridgen() const {
    soca_geo_gridgen_f90(keyGeom_);
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
  atlas::FunctionSpace * Geometry::atlasFunctionSpace(const std::string & functionSpaceName,
    bool halo) const {
    if (halo) {
      // Return function space with halo
      return atlasFunctionSpaceHalo_.get();
    } else {
      // Return function space without halo
      return atlasFunctionSpace_.get();
    }
  }
  // -----------------------------------------------------------------------------
  atlas::FieldSet * Geometry::atlasFieldSet(const std::string & functionSpaceName) const {
    return atlasFieldSet_.get();
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
