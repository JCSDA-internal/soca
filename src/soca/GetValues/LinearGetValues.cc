/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/GetValues/LinearGetValues.h"
#include "soca/GetValues/GetValuesFortran.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
LinearGetValues::LinearGetValues(const Geometry & geom, const ufo::Locations & locs) : locs_(locs),
  geom_(new Geometry(geom))
  {
  oops::Log::trace() << "LinearGetValues::LinearGetValues starting" << std::endl;
  soca_getvalues_create_f90(keyLinearGetValues_, geom.toFortran(), locs.toFortran());
  oops::Log::trace() << "LinearGetValues::LinearGetValues done" << std::endl;
}

LinearGetValues::~LinearGetValues() {
  soca_getvalues_delete_f90(keyLinearGetValues_);
  oops::Log::trace() << "LinearGetValues::~LinearGetValues done" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interpolate to obs locations
// -----------------------------------------------------------------------------
void LinearGetValues::setTrajectory(const State & state,
                                    const util::DateTime & t1, const util::DateTime & t2,
                                    ufo::GeoVaLs & geovals) {
  oops::Log::trace() << "LinearGetValues::setTrajectory starting" << std::endl;
  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  soca_getvalues_fill_geovals_f90(keyLinearGetValues_,
                                  geom_->toFortran(),
                                  state.toFortran(),
                                  &t1p, &t2p,
                                  locs_.toFortran(),
                                  geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::setTrajectory done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsTL(const Increment & incr,
                                    const util::DateTime & t1, const util::DateTime & t2,
                                    ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeovalsTL starting" << std::endl;
  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  soca_getvalues_fill_geovals_tl_f90(keyLinearGetValues_,
                                     geom_->toFortran(),
                                     incr.toFortran(),
                                     &t1p, &t2p,
                                     locs_.toFortran(),
                                     geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::fillGeovalsTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsAD(Increment & incr,
                                    const util::DateTime & t1, const util::DateTime & t2,
                                    const ufo::GeoVaLs & geovals) const {
  oops::Log::trace() << "LinearGetValues::fillGeovalsAD starting" << std::endl;
  const util::DateTime * t1p = &t1;
  const util::DateTime * t2p = &t2;
  soca_getvalues_fill_geovals_ad_f90(keyLinearGetValues_,
                                     geom_->toFortran(),
                                     incr.toFortran(),
                                     &t1p, &t2p,
                                     locs_.toFortran(),
                                     geovals.toFortran());
  oops::Log::trace() << "LinearGetValues::fillGeovalsAD done" << std::endl;
}

// -----------------------------------------------------------------------------
void LinearGetValues::print(std::ostream & os) const {
  os << "LinearGetValues" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace soca
