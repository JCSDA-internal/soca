/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "soca/State/State.h"
#include "soca/Increment/Increment.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Metadata.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace soca {


// -----------------------------------------------------------------------------

class UnstructuredInterpolator : public util::Printable,
                                 private util::ObjectCounter<UnstructuredInterpolator> {
 public:
  static const std::string classname() {return "soca::UnstructuredInterpolator";}

  UnstructuredInterpolator(const eckit::Configuration &,
                           const std::vector<double> &,
                           const std::vector<double> &,
                           const std::vector<double> &);

  void apply(const oops::Variables &, const atlas::FieldSet &, std::vector<double> &) const;
  void applyAD(const oops::Variables &, atlas::FieldSet &, const std::vector<double> &) const;

 private:
  void apply1lev(const std::string &, const atlas::array::ArrayView<double, 1> &,
                 std::vector<double>::iterator &) const;
  void applyLevs(const std::string &, const atlas::array::ArrayView<double, 2> &,
                 std::vector<double>::iterator &, const size_t &) const;
  void apply1levAD(const std::string &, atlas::array::ArrayView<double, 1> &,
                   std::vector<double>::const_iterator &) const;
  void applyLevsAD(const std::string &, atlas::array::ArrayView<double, 2> &,
                   std::vector<double>::const_iterator &, const size_t &) const;
  void print(std::ostream &) const override;

  std::string interp_method_;
  int nninterp_;
  size_t nout_;
  std::vector<std::vector<size_t>> interp_i_;
  std::vector<std::vector<double>> interp_w_;
};

// -----------------------------------------------------------------------------

UnstructuredInterpolator::UnstructuredInterpolator(const eckit::Configuration & config,
                                                   const std::vector<double> & lats_in,
                                                   const std::vector<double> & lons_in,
                                                   const std::vector<double> & latlon_out)
  : interp_method_(), nninterp_(0), nout_(0), interp_i_(), interp_w_()
{
  oops::Log::trace() << "UnstructuredInterpolator::UnstructuredInterpolator start" << std::endl;
  util::Timer timer("soca::UnstructuredInterpolator", "UnstructuredInterpolator");

  const atlas::Geometry earth(atlas::util::Earth::radius());
  const double close = 1.0e-10;

  // Create local input grid kd-tree
  // TODO(Travis) construct our multiple kdtrees once at a higher level... If I care
  const size_t npoints = lats_in.size();
  std::vector<size_t> indx(npoints);
  for (size_t jj = 0; jj < npoints; ++jj) indx[jj] = jj;

  atlas::util::IndexKDTree localTree(earth);
  localTree.build(lons_in, lats_in, indx);

  // Compute weights
  nninterp_ = config.getInt("nnearest", 4);
  nout_ = latlon_out.size() / 2;
  ASSERT(latlon_out.size() == 2 * nout_);
  interp_i_.resize(nout_, std::vector<size_t>(nninterp_));
  interp_w_.resize(nout_, std::vector<double>(nninterp_, 0.0));

  // This is a new option for this class, so isn't in any YAMLs yet!
  interp_method_ = config.getString("interpolation method", "barycentric");
  ASSERT(interp_method_ == "barycentric" || interp_method_ == "inverse distance");

  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    const double lat = latlon_out[2 * jloc];
    const double lon = latlon_out[2 * jloc + 1];
    atlas::PointLonLat obsloc(lon, lat);
    obsloc.normalise();

    const atlas::util::KDTree<size_t>::ValueList neighbours =
                                  localTree.closestPoints(obsloc, nninterp_);

    // Barycentric and inverse-distance interpolation both rely on indices, 1/distances
    size_t jj = 0;
    for (const atlas::util::KDTree<size_t>::Value & val : neighbours) {
      interp_i_[jloc][jj] = val.payload();
      interp_w_[jloc][jj] = 1.0 / val.distance();
      ++jj;
    }
    ASSERT(jj == nninterp_);

    if (interp_w_[jloc][0] > 1e10) {
      // Handle edge case where output point is close to one input point => interp_w very large
      // Atlas returns the neighbors in nearest-first order, so only need to check first element
      for (size_t jn = 0; jn < nninterp_; ++jn) interp_w_[jloc][jn] = 0.0;
      interp_w_[jloc][0] = 1.0;
    } else if (interp_method_ == "barycentric") {
      // Barycentric weights
      std::vector<double> bw(nninterp_);
      for (size_t j = 0; j < nninterp_; ++j) {
        double wprod = 1.0;
        const auto& p1 = neighbours[j].point();
        for (size_t k = 0; k < nninterp_; ++k) {
          if (k != j) {
            const auto& p2 = neighbours[k].point();
            const double dist = earth.distance(p1, p2);
            wprod *= std::max(dist, close);
          }
        }
        bw[j] = 1.0 / wprod;
      }
      // Interpolation weights from barycentric weights
      double wsum = 0.0;
      for (size_t j = 0; j < nninterp_; ++j) wsum += interp_w_[jloc][j] * bw[j];
      for (size_t j = 0; j < nninterp_; ++j) {
        interp_w_[jloc][j] *= bw[j] / wsum;
        ASSERT(interp_w_[jloc][j] >= 0.0 && interp_w_[jloc][j] <= 1.0);
      }
    } else {
      // Inverse-distance interpolation weights
      double wsum = 0.0;
      for (size_t j = 0; j < nninterp_; ++j) wsum += interp_w_[jloc][j];
      for (size_t j = 0; j < nninterp_; ++j) {
        interp_w_[jloc][j] /= wsum;
        ASSERT(interp_w_[jloc][j] >= 0.0 && interp_w_[jloc][j] <= 1.0);
      }
    }
  }
  oops::Log::trace() << "UnstructuredInterpolator::UnstructuredInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::apply(const oops::Variables & vars, const atlas::FieldSet & fset,
                                     std::vector<double> & vals) const {
  oops::Log::trace() << "UnstructuredInterpolator::apply starting" << std::endl;
  util::Timer timer("soca::UnstructuredInterpolator", "apply");

  size_t nflds = 0;
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string fname = vars[jf];
    atlas::Field fld = fset.field(fname);
    const size_t rank = fld.rank();
    ASSERT(rank >= 1 && rank <= 2);
    if (rank == 1) {
      nflds += 1;
    } else {
      nflds += fld.levels();
    }
  }
  vals.resize(nout_ * nflds);

  auto current = vals.begin();
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string fname = vars[jf];
    atlas::Field fld = fset.field(fname);
    const size_t rank = fld.rank();

    const std::string interp_type = fld.metadata().get<std::string>("interp_type");
    ASSERT(interp_type == "default" || interp_type == "integer" || interp_type == "nearest");

    if (rank == 1) {
      const atlas::array::ArrayView<double, 1> fldin = atlas::array::make_view<double, 1>(fld);
      this->apply1lev(interp_type, fldin, current);
    } else {
      const atlas::array::ArrayView<double, 2> fldin = atlas::array::make_view<double, 2>(fld);
      for (size_t jlev = 0; jlev < fldin.shape(1); ++jlev) {
        this->applyLevs(interp_type, fldin, current, jlev);
      }
    }
  }
  oops::Log::trace() << "UnstructuredInterpolator::apply done" << std::endl;
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::applyAD(const oops::Variables & vars, atlas::FieldSet &fset,
                                       const std::vector<double> & vals) const {
  oops::Log::trace() << "UnstructuredInterpolator::applyAD starting" << std::endl;
  util::Timer timer("soca::UnstructuredInterpolator", "applyAD");

  std::vector<double>::const_iterator current = vals.begin();
  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string fname = vars[jf];
    atlas::Field fld = fset.field(fname);
    const size_t rank = fld.rank();
    ASSERT(rank >= 1 && rank <= 2);

//    const std::string interp_type = fld.metadata().get<std::string>("interp_type");
//    ASSERT(interp_type == "default" || interp_type == "integer" || interp_type == "nearest");
    const std::string interp_type = "default";

    if (rank== 1) {
      atlas::array::ArrayView<double, 1> fldin = atlas::array::make_view<double, 1>(fld);
      this->apply1levAD(interp_type, fldin, current);
    } else {
      atlas::array::ArrayView<double, 2> fldin = atlas::array::make_view<double, 2>(fld);
      for (size_t jlev = 0; jlev < fldin.shape(1); ++jlev) {
        this->applyLevsAD(interp_type, fldin, current, jlev);
      }
    }
  }
  oops::Log::trace() << "UnstructuredInterpolator::applyAD done" << std::endl;
}


// -----------------------------------------------------------------------------

void UnstructuredInterpolator::apply1lev(const std::string & interp_type,
                                         const atlas::array::ArrayView<double, 1> & gridin,
                                         std::vector<double>::iterator & gridout) const {
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    *gridout = 0.0;
    if (interp_type == "default") {
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        *gridout += interp_w_[jloc][jj] * gridin(interp_i_[jloc][jj]);
      }
    } else if (interp_type == "integer") {
      // Find which integer value has largest weight in the stencil. We do this by taking two
      // passes through the (usually short) data: first to identify range of values, then to
      // determine weights for each integer.
      // Note that a std::map would be shorter to code, because it would avoid needing to find
      // the range of possible integer values, but vectors are almost always much more efficient.
      int minval = std::numeric_limits<int>().max();
      int maxval = std::numeric_limits<int>().min();
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        minval = std::min(minval, static_cast<int>(std::round(gridin(interp_i_[jloc][jj]))));
        maxval = std::max(maxval, static_cast<int>(std::round(gridin(interp_i_[jloc][jj]))));
      }
      std::vector<double> int_weights(maxval - minval + 1, 0.0);
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        const int this_int = std::round(gridin(interp_i_[jloc][jj]));
        int_weights[this_int - minval] += interp_w_[jloc][jj];
      }
      *gridout = minval + std::distance(int_weights.begin(),
          std::max_element(int_weights.begin(), int_weights.end()));
    } else if (interp_type == "nearest") {
      *gridout = gridin(interp_i_[jloc][0]);
    } else {
      throw eckit::BadValue("Unknown interpolation type");
    }
    ++gridout;
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::applyLevs(const std::string & interp_type,
                                         const atlas::array::ArrayView<double, 2> & gridin,
                                         std::vector<double>::iterator & gridout,
                                         const size_t & ilev) const {
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    *gridout = 0.0;
    if (interp_type == "default") {
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        *gridout += interp_w_[jloc][jj] * gridin(interp_i_[jloc][jj], ilev);
      }
    } else if (interp_type == "integer") {
      int minval = std::numeric_limits<int>().max();
      int maxval = std::numeric_limits<int>().min();
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        minval = std::min(minval, static_cast<int>(std::round(gridin(interp_i_[jloc][jj], ilev))));
        maxval = std::max(maxval, static_cast<int>(std::round(gridin(interp_i_[jloc][jj], ilev))));
      }
      std::vector<double> int_weights(maxval - minval + 1, 0.0);
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        const int this_int = std::round(gridin(interp_i_[jloc][jj], ilev));
        int_weights[this_int - minval] += interp_w_[jloc][jj];
      }
      *gridout = minval + std::distance(int_weights.begin(),
          std::max_element(int_weights.begin(), int_weights.end()));
    } else if (interp_type == "nearest") {
      *gridout = gridin(interp_i_[jloc][0], ilev);
    } else {
      throw eckit::BadValue("Unknown interpolation type");
    }
    ++gridout;
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::apply1levAD(const std::string & interp_type,
                                           atlas::array::ArrayView<double, 1> & gridin,
                                           std::vector<double>::const_iterator & gridout) const {
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    if (interp_type == "default") {
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        gridin(interp_i_[jloc][jj]) += interp_w_[jloc][jj] * *gridout;
      }
    } else if (interp_type == "integer") {
      throw eckit::BadValue("No adjoint for integer interpolation");
    } else if (interp_type == "nearest") {
      gridin(interp_i_[jloc][0]) += *gridout;
    } else {
      throw eckit::BadValue("Unknown interpolation type");
    }
    ++gridout;
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::applyLevsAD(const std::string & interp_type,
                                           atlas::array::ArrayView<double, 2> & gridin,
                                           std::vector<double>::const_iterator & gridout,
                                           const size_t & ilev) const {
  for (size_t jloc = 0; jloc < nout_; ++jloc) {
    if (interp_type == "default") {
      for (size_t jj = 0; jj < nninterp_; ++jj) {
        gridin(interp_i_[jloc][jj], ilev) += interp_w_[jloc][jj] * *gridout;
      }
    } else if (interp_type == "integer") {
      throw eckit::BadValue("No adjoint for integer interpolation");
    } else if (interp_type == "nearest") {
      gridin(interp_i_[jloc][0], ilev) += *gridout;
    } else {
      throw eckit::BadValue("Unknown interpolation type");
    }
    ++gridout;
  }
}

// -----------------------------------------------------------------------------

void UnstructuredInterpolator::print(std::ostream & os) const
{
  os << "SOCAUnstructuredInterpolator";
}

// -----------------------------------------------------------------------------

}  // namespace soca
