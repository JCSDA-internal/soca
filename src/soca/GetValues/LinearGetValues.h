/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GETVALUES_LINEARGETVALUES_H_
#define SOCA_GETVALUES_LINEARGETVALUES_H_

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "soca/Fortran.h"
#include "soca/Increment/Increment.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "ufo/Locations.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace ufo {
  class GeoVaLs;
  class Locations;
}
namespace soca {
  class Geometry;
  class State;
}

//-----------------------------------------------------------------------------

namespace soca {

  /// SOCA GetValues
  /*!
   * GetValues interpolates State to observation locations
   */
class LinearGetValues : public util::Printable,
                        private util::ObjectCounter<LinearGetValues> {
 public:
  static const std::string classname() {return "soca::LinearGetValues";}

  ///  Constructor, destructor
  LinearGetValues(const Geometry &, const ufo::Locations &);
  virtual ~LinearGetValues();

  /// Trajectory for the linearized interpolation
  void setTrajectory(const State & state, const util::DateTime & t1, const util::DateTime & t2,
                     ufo::GeoVaLs & geovals);

  /// Forward and backward interpolation
  void fillGeoVaLsTL(const Increment & inc,
                     const util::DateTime & t1, const util::DateTime & t2,
                     ufo::GeoVaLs & geovals) const;
  void fillGeoVaLsAD(Increment & inc,
                     const util::DateTime & t1, const util::DateTime & t2,
                     const ufo::GeoVaLs & geovals) const;

 private:
  void print(std::ostream &) const;
  F90getval keyLinearGetValues_;
  ufo::Locations locs_;
  std::shared_ptr<const Geometry> geom_;
};
// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_GETVALUES_LINEARGETVALUES_H_
