/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/DateTime.h"

#include "soca/Geometry/Geometry.h"
#include "soca/GetValues/GetValuesFortran.h"
#include "soca/GetValues/LinearGetValues.h"
#include "soca/Increment/Increment.h"
#include "soca/State/State.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace soca {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
LinearGetValues::LinearGetValues(const Geometry & geom,
                                 const ufo::Locations & locs)
  : locs_(locs), geom_(new Geometry(geom))
{
  soca_getvalues_create_f90(keyLinearGetValues_,
                            geom.toFortran(),
                            locs.toFortran());
}

LinearGetValues::~LinearGetValues()
{
  soca_getvalues_delete_f90(keyLinearGetValues_);
}
// -----------------------------------------------------------------------------
/// Interpolate to obs locations
// -----------------------------------------------------------------------------
void LinearGetValues::setTrajectory(const State & state,
                                    const util::DateTime & t1,
                                    const util::DateTime & t2,
                                    ufo::GeoVaLs & geovals) {
  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  soca_getvalues_fill_geovals_f90(keyLinearGetValues_,
                                  geom_->toFortran(),
                                  state.toFortran(),
                                  &t1p, &t2p,
                                  locs_.toFortran(),
                                  geovals.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsTL(const Increment & incr,
                                    const util::DateTime & t1,
                                    const util::DateTime & t2,
                                    ufo::GeoVaLs & geovals) const {
  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  soca_getvalues_fill_geovals_tl_f90(keyLinearGetValues_,
                                     geom_->toFortran(),
                                     incr.toFortran(),
                                     &t1p, &t2p,
                                     locs_.toFortran(),
                                     geovals.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsAD(Increment & incr,
                                    const util::DateTime & t1,
                                    const util::DateTime & t2,
                                    const ufo::GeoVaLs & geovals) const {
  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  soca_getvalues_fill_geovals_ad_f90(keyLinearGetValues_,
                                     geom_->toFortran(),
                                     incr.toFortran(),
                                     &t1p, &t2p,
                                     locs_.toFortran(),
                                     geovals.toFortran());
}

// -----------------------------------------------------------------------------
void LinearGetValues::print(std::ostream & os) const {
  os << "LinearGetValues" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace soca
