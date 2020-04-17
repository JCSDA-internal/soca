/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GETVALUES_GETVALUES_H_
#define SOCA_GETVALUES_GETVALUES_H_

#include <ostream>
#include <string>
#include <memory>

#include "soca/Fortran.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/Locations.h"

// Forward declarations
namespace ufo {
  class GeoVaLs;
}
namespace soca {
  class Geometry;
  class State;
}

//-----------------------------------------------------------------------------

namespace soca {

  /// SOCA GetValues
  /*!
   * GetValues: interpolate State to observation locations
   */
class GetValues : public util::Printable,
                    private util::ObjectCounter<GetValues> {
 public:
  static const std::string classname() {return "soca::GetValues";}

/// saves all locations locs to use during filling GeoVaLs
  GetValues(const Geometry &, const ufo::Locations & locs);
  virtual ~GetValues();

  /// fills in geovals for all observations in the timeframe (t1, t2],
  /// geovals are interpolated trilinearly from state at the nearest gridpoints
  void fillGeoVaLs(const State &,
                   const util::DateTime & t1,
                   const util::DateTime & t2,
                   ufo::GeoVaLs &) const;

  /// Read interpolated GeoVaLs at observation location
  void getValuesFromFile(const ufo::Locations &,
                         const oops::Variables &,
                         ufo::GeoVaLs &) const;

 private:
  void print(std::ostream &) const;
  F90getval keyGetValues_;
  ufo::Locations locs_;
  std::shared_ptr<const Geometry> geom_;
};
// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_GETVALUES_GETVALUES_H_
