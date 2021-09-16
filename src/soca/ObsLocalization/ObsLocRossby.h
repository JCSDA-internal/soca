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

#include "ufo/obslocalization/ObsLocalization.h"
#include "ufo/obslocalization/ObsLocParameters.h"

#include "soca/ObsLocalization/ObsLocRossbyParameters.h"

// -----------------------------------------------------------------------------

namespace soca {

template<class MODEL>
class ObsLocRossby: public ufo::ObsLocalization<MODEL> {
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
    ufo::ObsLocalization<MODEL>::ObsLocalization(config, obsspace) {
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

  // do distance search
  ufo::ObsLocalization<MODEL>::searchObs(i, lengthscale);

  // return refs to internals of ObsLocalization
  const std::vector<int> & localobs = ufo::ObsLocalization<MODEL>::localobs();
  const std::vector<double> & horizontalObsdist =
    ufo::ObsLocalization<MODEL>::horizontalObsdist();

  // set all to missing (outside of localization distance)
  const double missing = util::missingValue(double());
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    locvector[jj] = missing;
  }

  // calculate distance based localization
  const size_t nvars = locvector.nvars();
  for (size_t jlocal = 0; jlocal < localobs.size(); ++jlocal) {
    double locFactor = oops::gc99(horizontalObsdist[jlocal] / lengthscale);
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      locvector[jvar + localobs[jlocal] * nvars] = locFactor;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
