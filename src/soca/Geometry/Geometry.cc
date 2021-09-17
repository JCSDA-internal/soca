/*
 * (C) Copyright 2017-2020 UCAR
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
      atmconf_(conf),
      initatm_(initAtm(conf)),
      fmsinput_(comm, conf) {

    fmsinput_.updateNameList();

    soca_geo_setup_f90(keyGeom_, &conf, &comm);

    // Set ATLAS lon/lat field
    atlasFieldSet_.reset(new atlas::FieldSet());
    soca_geo_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet_->get());
    atlas::Field atlasField = atlasFieldSet_->field("lonlat");

    // Create ATLAS function space
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));

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
      atmconf_(other.atmconf_),
      initatm_(initAtm(other.atmconf_)),
      fmsinput_(other.fmsinput_) {
    const int key_geo = other.keyGeom_;
    soca_geo_clone_f90(keyGeom_, key_geo);
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                              other.atlasFunctionSpace_->lonlat()));
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
    int ist, iend, jst, jend, kst, kend;
    soca_geo_start_end_f90(keyGeom_, ist, iend, jst, jend, kst, kend);
    return GeometryIterator(*this, ist, jst, kst);
  }
  // -----------------------------------------------------------------------------
  GeometryIterator Geometry::end() const {
    // return end of the geometry on this mpi tile
    // decided to return index out of bounds for the iterator loops to work
    return GeometryIterator(*this, -1, -1, -1);
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
  atlas::FunctionSpace * Geometry::atlasFunctionSpace() const {
    return atlasFunctionSpace_.get();
  }
  // -----------------------------------------------------------------------------
  atlas::FieldSet * Geometry::atlasFieldSet() const {
    return atlasFieldSet_.get();
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
