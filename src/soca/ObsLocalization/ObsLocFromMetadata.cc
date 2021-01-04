/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/ObsLocalization/ObsLocFromMetadata.h"
#include "oops/generic/gc99.h"
#include "oops/interface/ObsLocalization.h"
#include "ufo/ObsTraits.h"

namespace soca {

// -----------------------------------------------------------------------------

oops::ObsLocalizationMaker<
       ufo::ObsTraits,
       oops::ObsLocalization<ufo::ObsTraits, ObsLocFromMetadata> >
  ObsLocFromMetadata_("obs metadata");

// -----------------------------------------------------------------------------

ObsLocFromMetadata::ObsLocFromMetadata(const eckit::Configuration & config,
                           const ioda::ObsSpace & obsdb)
  : obsdb_(obsdb) {}

// -----------------------------------------------------------------------------

void ObsLocFromMetadata::multiply(ioda::ObsVector & dy) const {
  const std::vector<double> & obsdist = obsdb_.obsdist();
  const size_t nlocs = dy.nlocs();
  const size_t nvars = dy.nvars();

  std::vector<double> loc(nlocs);
  obsdb_.get_db("MetaData", "obs_localization", loc);

  // create a temporary vector with localization values
  // (by putting in an ObsVector, multiplication respects the missing values)
  ioda::ObsVector gcVec(dy);
  gcVec.zero();
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    double gc = oops::gc99(obsdist[jloc] / (loc[jloc] * 2.0 / sqrt(0.3)));
    for (size_t jvar = 0; jvar < nvars; ++jvar) {
      gcVec[jvar + jloc * nvars] = gc;
    }
  }

  dy *= gcVec;
}

// -----------------------------------------------------------------------------

}  // namespace soca
