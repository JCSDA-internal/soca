/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
#define SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_

#include <algorithm>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/generic/gc99.h"

#include "ufo/obslocalization/ObsLocGC99.h"
#include "ufo/obslocalization/ObsLocParameters.h"

#include "soca/ObsLocalization/ObsLocRossbyParameters.h"

// -----------------------------------------------------------------------------

namespace soca {

template<class MODEL>
class ObsLocRossby: public ufo::ObsLocGC99<MODEL> {
  typedef typename MODEL::GeometryIterator   GeometryIterator_;

 public:
  ObsLocRossby(const eckit::Configuration &, const ioda::ObsSpace &);
  void computeLocalization(
    const GeometryIterator_ &,
    ioda::ObsVector & locfactor) const override;

 private:
  ObsLocRossbyParameters options_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
ObsLocRossby<MODEL>::ObsLocRossby(
      const eckit::Configuration& config,
      const ioda::ObsSpace & obsspace):
    ufo::ObsLocGC99<MODEL>::ObsLocGC99(config, obsspace) {
  options_.deserialize(config);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsLocRossby<MODEL>::computeLocalization(
    const GeometryIterator_ & i,
    ioda::ObsVector & locvector) const {

  // calculate the length scale at this location
  double lengthscale = options_.base;
  lengthscale += options_.mult * i.getRossbyRadius();
  lengthscale = std::max(lengthscale, options_.min_grid * sqrt(i.getArea()));
  const boost::optional<double> & minval = options_.min;
  const boost::optional<double> & maxval = options_.max;
  if (minval != boost::none) lengthscale = std::max(lengthscale, *minval);
  if (maxval != boost::none) lengthscale = std::min(lengthscale, *maxval);

  // convert from gaussian to gaspari-cohn width
  lengthscale *= 2.0/sqrt(0.3);

  // do GC99 localization
  ufo::ObsLocGC99<MODEL>::setLengthscale(lengthscale);
  ufo::ObsLocGC99<MODEL>::computeLocalization(i, locvector);
}

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
