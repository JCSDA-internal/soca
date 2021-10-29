/*
 * (C) Copyright 2020-2021 UCAR
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
#include "soca/Transforms/Model2GeoVaLs/LinearModel2GeoVaLs.h"
#include "soca/Transforms/Model2GeoVaLs/Model2GeoVaLs.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace soca {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
LinearGetValues::LinearGetValues(const Geometry & geom,
                                 const ufo::Locations & locs,
                                 const eckit::Configuration &config)
  : locs_(locs), geom_(new Geometry(geom)),
    model2geovals_(new Model2GeoVaLs(geom, config))
{
  soca_getvalues_create_f90(keyLinearGetValues_,
                            geom.toFortran(),
                            locs);
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
  std::unique_ptr<State> varChangeState;
  const State * state_ptr;

  // Do variable change if it has not already been done.
  // TODO(travis): remove this once Yannick is done rearranging things in oops.
  if ( geovals.getVars() <= state.variables() ) {
    state_ptr = &state;
  } else {
    varChangeState.reset(new State(*geom_, geovals.getVars(),
                         state.validTime()));
    model2geovals_->changeVar(state, *varChangeState);
    state_ptr = varChangeState.get();
  }

  // TODO(travis) : change to a map to store multiple time slices?
  eckit::LocalConfiguration conf;
  linearmodel2geovals_.reset(new LinearModel2GeoVaLs(state, state,
                                                     *geom_, conf));

  soca_getvalues_fill_geovals_f90(keyLinearGetValues_,
                                  geom_->toFortran(),
                                  state_ptr->toFortran(),
                                  t1, t2, locs_,
                                  geovals.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsTL(const Increment & incr,
                                    const util::DateTime & t1,
                                    const util::DateTime & t2,
                                    ufo::GeoVaLs & geovals) const {
  Increment incrGeovals(*geom_, geovals.getVars(), incr.validTime());
  linearmodel2geovals_->multiply(incr, incrGeovals);
  soca_getvalues_fill_geovals_tl_f90(keyLinearGetValues_,
                                     geom_->toFortran(),
                                     incrGeovals.toFortran(),
                                     t1, t2, locs_,
                                     geovals.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinearGetValues::fillGeoVaLsAD(Increment & incr,
                                    const util::DateTime & t1,
                                    const util::DateTime & t2,
                                    const ufo::GeoVaLs & geovals) const {
  Increment incrGeovals(*geom_, geovals.getVars(), incr.validTime());
  soca_getvalues_fill_geovals_ad_f90(keyLinearGetValues_,
                                     geom_->toFortran(),
                                     incrGeovals.toFortran(),
                                     t1, t2, locs_,
                                     geovals.toFortran());
  linearmodel2geovals_->multiplyAD(incrGeovals, incr);
}

// -----------------------------------------------------------------------------
void LinearGetValues::print(std::ostream & os) const {
  os << "LinearGetValues" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace soca
