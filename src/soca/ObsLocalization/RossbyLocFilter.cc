/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "soca/ObsLocalization/RossbyLocFilter.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/ObsTraits.h"

namespace soca {

// -----------------------------------------------------------------------------

oops::FilterMaker<
      ufo::ObsTraits,
      oops::ObsFilter<ufo::ObsTraits, RossbyLocFilter> >
  RossbyLocFilterMaker_("Rossby Localization");

// -----------------------------------------------------------------------------

RossbyLocFilter::RossbyLocFilter(ioda::ObsSpace & obsdb,
                           const Parameters_ & params,
                           std::shared_ptr<ioda::ObsDataVector<int> > flags,
                           std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, params, flags, obserr), parameters_(params)
{
  allvars_+= ufo::Variable("rossby_radius@GeoVaLs");
}

void RossbyLocFilter::applyFilter(const std::vector<bool> & apply,
                                const ufo::Variables & filtervars,
                                std::vector<std::vector<bool>> & flagged) const
{
  const size_t nlocs = obsdb_.nlocs();
  std::vector<float> loc(nlocs);
  loc.assign(nlocs, 1.0);

  data_.get(ufo::Variable("rossby_radius@GeoVaLs"), loc);
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs]) {
      loc[jobs] *= parameters_.multiplier.value();
      loc[jobs] = std::max(loc[jobs], parameters_.minvalue.value());
      loc[jobs] = std::min(loc[jobs], parameters_.maxvalue.value());
    }
  }
  obsdb_.put_db("MetaData", "obs_localization", loc);
}

}  // namespace soca
