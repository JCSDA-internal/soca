/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "eckit/config/YAMLConfiguration.h"

#include "soca/Geometry/Geometry.h"

// -----------------------------------------------------------------------------
namespace soca {

  const std::vector<char> grids{ 'h', 'u', 'v'};

  // -----------------------------------------------------------------------------
  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm),
      fmsinput_(comm, conf) {

    fmsinput_.updateNameList();

    soca_geo_setup_f90(keyGeom_, &conf, &comm);

    // Set ATLAS lonlat and function space (with and without halos)
    atlas::FieldSet lonlat;
    soca_geo_lonlat_f90(keyGeom_, lonlat.get());
    functionSpace_ = atlas::functionspace::PointCloud(lonlat->field("lonlat"));
    functionSpaceIncHalo_ = atlas::functionspace::PointCloud(lonlat->field("lonlat_inc_halos"));

    // Set ATLAS function space pointer in Fortran
    soca_geo_set_atlas_functionspace_pointer_f90(
      keyGeom_, functionSpace_.get(), functionSpaceIncHalo_.get());

    // Fill ATLAS fieldset
    soca_geo_to_fieldset_f90(keyGeom_, extraFields_.get());

    // create kdtrees
    int kdidx = 0;
    for (auto grid : grids) {
      std::vector<double> lats;
      std::vector<double> lons;
      size_t npoints;
      std::vector<size_t> indx;

      this->latlon(lats, lons, true, grid, false);
      npoints = lats.size();
      indx.resize(npoints);
      for (size_t jj = 0; jj < npoints; ++jj) indx[jj] = jj;
      localTree_[kdidx++].build(lons, lats, indx);

      this->latlon(lats, lons, true, grid, true);
      npoints = lats.size();
      indx.resize(npoints);
      for (size_t jj = 0; jj < npoints; ++jj) indx[jj] = jj;
      localTree_[kdidx++].build(lons, lats, indx);
    }
  }
  // -----------------------------------------------------------------------------
  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_),
      fmsinput_(other.fmsinput_)
       {
    const int key_geo = other.keyGeom_;
    soca_geo_clone_f90(keyGeom_, key_geo);

    functionSpace_ = atlas::functionspace::PointCloud(other.functionSpace_->lonlat());
    functionSpaceIncHalo_ = atlas::functionspace::PointCloud(other.functionSpaceIncHalo_->lonlat());
    soca_geo_set_atlas_functionspace_pointer_f90(keyGeom_, functionSpace_.get(),
                                                 functionSpaceIncHalo_.get());

    extraFields_ = atlas::FieldSet();
    for (int jfield = 0; jfield < other.extraFields_->size(); ++jfield) {
      atlas::Field atlasField = other.extraFields_->field(jfield);
      extraFields_->add(atlasField);
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
  void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
      const bool halo) const {
    // Assume that we can get by with just using the unmasked H grid here (i.e.
    // ignore the different staggered grids)
    // (This method should only be called by code in OOPS used for determining the
    //  PE that owns a specific point.... so it doesn't have to be exact)
    latlon(lats, lons, halo, 'h', false);
  }
  // -----------------------------------------------------------------------------
  void Geometry::latlon(
        std::vector<double> & lats, std::vector<double> & lons, const bool halo,
        const char grid, const bool masked) const {
    // get the number of gridpoints
    int gridSize;
    soca_geo_gridsize_f90(keyGeom_, grid, masked, halo, gridSize);

    // get the lat/lon of those gridpoints
    lats.resize(gridSize);
    lons.resize(gridSize);
    soca_geo_gridlatlon_f90(keyGeom_, grid, masked, halo, gridSize,
      lats.data(), lons.data());
  }
  // -----------------------------------------------------------------------------
  atlas::util::KDTree<size_t>::ValueList Geometry::closestPoints(
    const double lat, const double lon, const int npoints,
    const char grid, const bool masked) const {

    // determine which kdtree to use
    // TODO(travis) make sure not -1
    int kdidx = -1;
    for (int j=0; j < grids.size(); j++) {
      if (grids[j] == grid) kdidx = j*2;
    }
    if (masked) kdidx += 1;

    atlas::PointLonLat obsloc(lon, lat);
    obsloc.normalise();
    atlas::util::KDTree<size_t>::ValueList neighbours =
      localTree_[kdidx].closestPoints(obsloc, npoints);

    return neighbours;
    }
  // -----------------------------------------------------------------------------
  // Get the properties of the grid for a given variable (masked and u/v/h grid)
  void Geometry::getVarGrid(const std::string &var, char & grid, bool & masked) const {
      oops::Variables vars;
      vars.push_back(var);

      soca_geo_getvargrid_f90(keyGeom_, vars, grid, masked);
  }

}  // namespace soca
