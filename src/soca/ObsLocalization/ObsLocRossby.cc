/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/ObsLocalization/ObsLocRossby.h"

#include <algorithm>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/ObsLocalization/ObsLocRossbyParameters.h"
#include "ufo/obslocalization/ObsHorLocGC99.h"
#include "ufo/obslocalization/ObsLocalizationBase.h"

namespace soca {

// -----------------------------------------------------------------------------
static ufo::ObsLocalizationMaker<GeometryIterator, soca::ObsLocRossby> obslocrossby_("Rossby");
// -----------------------------------------------------------------------------

ObsLocRossby::ObsLocRossby(const eckit::Configuration & config,
                           ioda::ObsSpace & obsspace):
    ufo::ObsHorLocGC99<GeometryIterator>::ObsHorLocGC99(config, obsspace),
    options_(), cacheVector_(obsspace), cachePoint_(-999, -999) {
    options_.validateAndDeserialize(config);
}

// -----------------------------------------------------------------------------

void ObsLocRossby::computeLocalization(
    const GeometryIterator & i,
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
  ufo::ObsHorLocGC99<GeometryIterator>::getLocalObs(i, lengthscale);
  ufo::ObsHorLocGC99<GeometryIterator>::localizeLocalObs(i, locvector, localobs);
  cacheVector_ = locvector;
  cachePoint_ = refPoint2;
}

// -----------------------------------------------------------------------------

}  // namespace soca
