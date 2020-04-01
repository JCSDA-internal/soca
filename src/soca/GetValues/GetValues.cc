/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/GetValues/GetValues.h"
#include "soca/GetValues/GetValuesFortran.h"
#include "soca/State/State.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
GetValues::GetValues(const Geometry & geom, const ufo::Locations & locs) :
    locs_(locs), geom_(new Geometry(geom)) {
  //soca_getvalues_create_f90(key_, geom.toFortran(), locs.toFortran());
}
// -----------------------------------------------------------------------------
GetValues::~GetValues() {
  //soca_getvalues_delete_f90(key_);
}
// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void GetValues::fillGeoVaLs(const State & state,
                            const util::DateTime & t1, const util::DateTime & t2,
                            ufo::GeoVaLs & gom) const {
  const util::DateTime * pt1 = &t1;
  const util::DateTime * pt2 = &t2;
  //soca_getvalues_interp_f90(key_, state.fields().toFortran(),
  //                        &pt1, &pt2, gom.toFortran());
}
// -----------------------------------------------------------------------------
void GetValues::print(std::ostream & os) const {
  os << "GetValues" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace soca
