/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/ObsLocalization/ObsLocRossby.h"
#include "oops/generic/gc99.h"
#include "oops/interface/ObsLocalization.h"
#include "ufo/ObsTraits.h"

namespace soca {

// -----------------------------------------------------------------------------

oops::ObsLocalizationMaker<
       ufo::ObsTraits,
       oops::ObsLocalization<ufo::ObsTraits, ObsLocRossby> >
  maker_("rossby radius");

// -----------------------------------------------------------------------------

ObsLocRossby::ObsLocRossby(const eckit::Configuration & config,
                           const ioda::ObsSpace & obsdb)
  : obsdb_(obsdb),
    rscale_(config.getDouble("lengthscale")),
    locMult_(config.getDouble("loc mult")),
    locMin_(config.getDouble("loc min"))
  {}

// -----------------------------------------------------------------------------

void ObsLocRossby::multiply(ioda::ObsVector & dy) const {
  const std::vector<double> & obsdist = obsdb_.obsdist();
  const size_t nlocs = dy.nlocs();
  const size_t nvars = dy.nvars();

  // create a temporary vector with localization values
  // (by putting in an ObsVector, multiplication respects the missing values)
  ioda::ObsVector gcVec(dy);
  gcVec.zero();
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    double gc = oops::gc99(obsdist[jloc] / rscale_);
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      gcVec[jvar + jloc * nvars] = gc;
    }
  }

  dy *= gcVec;
}

// -----------------------------------------------------------------------------

}  // namespace soca
