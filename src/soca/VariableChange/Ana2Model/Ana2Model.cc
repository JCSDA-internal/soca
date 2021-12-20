/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"
#include "soca/VariableChange/Ana2Model/Ana2Model.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------

static VariableChangeMaker<Ana2Model> makerVarChaA2M_("Ana2Model");

// -----------------------------------------------------------------------------
Ana2Model::Ana2Model(const Geometry & resol, const eckit::Configuration & conf)
: uvars_(initRotate(conf, "u")), vvars_(initRotate(conf, "v")),
  logvars_(initTrans(conf, "var"))
{
  Log::trace() << "Ana2Model::Ana2Model start" << std::endl;
  ASSERT(uvars_.size() == vvars_.size());
  Log::trace() << "Ana2Model::Ana2Model Rotating:"
               << " u = " << uvars_ << " v = " << vvars_ << std::endl;
  Log::trace() << "Ana2Model::Log Transforming:"
               << " var = " << logvars_ << std::endl;
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
  xm.logtrans(logvars_);
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
  xa.expontrans(logvars_);
  xa.validTime() = xm.validTime();
  Log::trace() << "Ana2Model::changeVarInverse done" << xa << std::endl;
}
// -----------------------------------------------------------------------------
void Ana2Model::print(std::ostream & os) const {
  os << "Ana2Model";
}
// -----------------------------------------------------------------------------
}  // namespace soca
