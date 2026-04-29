/*
 * (C) Copyright 2021-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
#define SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

#include "ufo/obslocalization/ObsHorLocGC99.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/ObsLocalization/ObsLocRossbyParameters.h"

// -----------------------------------------------------------------------------

namespace soca {

class ObsLocRossby: public ufo::ObsHorLocGC99<GeometryIterator> {
  typedef typename ufo::ObsHorLocalization<GeometryIterator>::LocalObs LocalObs_;

 public:
  ObsLocRossby(const eckit::Configuration &, ioda::ObsSpace &);

  /// Compute Rossby radius based localization and update localization values
  /// in \p locvector. Missing values indicate that observation is outside of
  /// localization.
  void computeLocalization(
    const GeometryIterator &,
    ioda::ObsVector & locfactor) const override;

  double computeLocalization(const eckit::geometry::Point3 &,
                             const eckit::geometry::Point3 &) const override;

 private:
  ObsLocRossbyParameters options_;
  mutable ioda::ObsVector cacheVector_;
  mutable eckit::geometry::Point2 cachePoint_;
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_OBSLOCALIZATION_OBSLOCROSSBY_H_
