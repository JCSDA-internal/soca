/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsSpace.h"

#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"

#include "soca/Geometry/Geometry.h"
#include "soca/GetValues/GetValues.h"
#include "soca/GetValues/GetValuesFortran.h"
#include "soca/State/State.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace soca {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
GetValues::GetValues(const Geometry & geom,
                     const ufo::Locations & locs,
                     const eckit::Configuration & config)
  : locs_(locs), geom_(new Geometry(geom)) {
  soca_getvalues_create_f90(keyGetValues_, geom.toFortran(), locs);
}
// -----------------------------------------------------------------------------
GetValues::~GetValues()
{
  soca_getvalues_delete_f90(keyGetValues_);
}
// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void GetValues::fillGeoVaLs(const State & state,
                            const util::DateTime & t1,
                            const util::DateTime & t2,
                            ufo::GeoVaLs & geovals) const {
  soca_getvalues_fill_geovals_f90(keyGetValues_,
                                  geom_->toFortran(),
                                  state.toFortran(),
                                  t1, t2, locs_,
                                  geovals.toFortran());
}
// -----------------------------------------------------------------------------
void GetValues::print(std::ostream & os) const {
  os << "GetValues" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace soca
