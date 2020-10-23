/*
 * (C) Copyright 2020-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Transforms/LogExpon/LogExpon.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"
#include "oops/util/Logger.h"
#include "eckit/exception/Exceptions.h"

using oops::Log;

namespace soca {
// -----------------------------------------------------------------------------
LogExpon::LogExpon(const Geometry & resol, const eckit::Configuration & conf)
: trvars_(initTrans(conf, "var"))
{
  Log::trace() << "LogExpon::LogExpon start" << std::endl;
  Log::trace() << "LogExpon::LogExpon Transforming:"
               << " var = " << trvars_ << std::endl;
  Log::trace() << "LogExpon::LogExpon done" << std::endl;
}
// -----------------------------------------------------------------------------
LogExpon::~LogExpon() {
  oops::Log::trace() << "ChangeSOCA destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void LogExpon::changeVar(const State & xa, State & xm) const {
  oops::Log::trace() << "LogExpon::changeVar starting" << xa <<
                        std::endl;
  util::DateTime * vtime = &xm.validTime();
  xm = xa;
  xm.logtrans(trvars_);
  xm.validTime() = xa.validTime();
  Log::trace() << "LogExpon::changeVar done" << xm << std::endl;
}
// -----------------------------------------------------------------------------
void LogExpon::changeVarInverse(const State & xm, State & xa) const {
  Log::trace() << "LogExpon::changeVarInverse starting" <<xm <<
                        std::endl;
  util::DateTime * vtime = &xa.validTime();
  xa = xm;
  xa.expontrans(trvars_);
  xa.validTime() = xm.validTime();
  Log::trace() << "LogExpon::changeVarInverse done" << xa << std::endl;
}
// -----------------------------------------------------------------------------
void LogExpon::print(std::ostream & os) const {
  os << "LogExpon";
}
// -----------------------------------------------------------------------------
}  // namespace soca
