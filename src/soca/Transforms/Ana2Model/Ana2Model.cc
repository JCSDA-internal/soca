/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Transforms/Ana2Model/Ana2Model.h"

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
Ana2Model::Ana2Model(const Geometry & resol, const eckit::Configuration & conf)
: uvars_(initRotate(conf, "u")), vvars_(initRotate(conf, "v"))
{
  Log::trace() << "Ana2Model::Ana2Model start" << std::endl;
  ASSERT(uvars_.size() == vvars_.size());
  Log::trace() << "Ana2Model::Ana2Model Rotating:"
               << " u = " << uvars_ << " v = " << vvars_ << std::endl;
  Log::trace() << "Ana2Model::Ana2Model done" << std::endl;
}
// -----------------------------------------------------------------------------
Ana2Model::~Ana2Model() {
  oops::Log::trace() << "ChangeSOCA destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void Ana2Model::changeVar(const State & xa,
                                State & xm) const {
  oops::Log::trace() << "Ana2Model::changeVar starting" << xa <<
                        std::endl;
  util::DateTime * vtime = &xm.validTime();
  xm = xa;
  xm.rotate2grid(uvars_, vvars_);
  xm.validTime() = xa.validTime();
  Log::trace() << "Ana2Model::changeVar done" << xm << std::endl;
}
// -----------------------------------------------------------------------------
void Ana2Model::changeVarInverse(const State & xm,
                                       State & xa) const {
  Log::trace() << "Ana2Model::changeVarInverse starting" <<xm <<
                        std::endl;
  util::DateTime * vtime = &xa.validTime();
  xa = xm;
  xa.rotate2north(uvars_, vvars_);
  xa.validTime() = xm.validTime();
  Log::trace() << "Ana2Model::changeVarInverse done" << xa << std::endl;
}
// -----------------------------------------------------------------------------
void Ana2Model::print(std::ostream & os) const {
  os << "Ana2Model";
}
// -----------------------------------------------------------------------------
}  // namespace soca
