/*
 * (C) Copyright 2021-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
#define SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_

#include <algorithm>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/obslocalization/ObsHorLocGC99.h"
#include "soca/ObsLocalization/ObsLocRossbyParameters.h"

// -----------------------------------------------------------------------------

namespace soca {

template<class MODEL>
class ObsLocRossby: public ufo::ObsHorLocGC99<MODEL> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;
  typedef typename ufo::ObsHorLocalization<MODEL>::LocalObs LocalObs_;

 public:
  ObsLocRossby(const eckit::Configuration &, ioda::ObsSpace &);

  /// Compute Rossby radius based localization and update localization values
  /// in \p locvector. Missing values indicate that observation is outside of
  /// localization.
  void computeLocalization(
    const GeometryIterator_ &,
    ioda::ObsVector & locfactor) const override;

 private:
  ObsLocRossbyParameters options_;
  mutable ioda::ObsVector cacheVector_;
  mutable eckit::geometry::Point2 cachePoint_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsLocRossby<MODEL>::ObsLocRossby(const eckit::Configuration & config,
                                  ioda::ObsSpace & obsspace):
    ufo::ObsHorLocGC99<MODEL>::ObsHorLocGC99(config, obsspace),
    options_(), cacheVector_(obsspace), cachePoint_(-999, -999) {
    options_.validateAndDeserialize(config);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocRossby<MODEL>::computeLocalization(
    const GeometryIterator_ & i,
    ioda::ObsVector & locvector) const {
  const eckit::geometry::Point3 refPoint = *i;
  const eckit::geometry::Point2 refPoint2(refPoint[0], refPoint[1]);
  if (refPoint2 == cachePoint_) {
    locvector = cacheVector_;
    return;
  }

  // calculate the length scale at this location
  double lengthscale = options_.base;
  lengthscale += options_.mult * i.getFieldValue("rossby_radius");
  lengthscale = std::max(lengthscale, options_.min_grid * sqrt(i.getFieldValue("area")));
  const boost::optional<double> & minval = options_.min;
  const boost::optional<double> & maxval = options_.max;
  if (minval != boost::none) lengthscale = std::max(lengthscale, *minval);
  if (maxval != boost::none) lengthscale = std::min(lengthscale, *maxval);

  // convert from gaussian to gaspari-cohn width
  lengthscale *= 2.0/sqrt(0.3);

  // Apply GC99 localization
  const LocalObs_ & localobs =
  ufo::ObsHorLocGC99<MODEL>::getLocalObs(i, lengthscale);
  ufo::ObsHorLocGC99<MODEL>::localizeLocalObs(i, locvector, localobs);
  cacheVector_ = locvector;
  cachePoint_ = refPoint2;
}

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
