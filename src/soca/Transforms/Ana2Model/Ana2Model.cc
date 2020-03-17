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

using oops::Log;

namespace soca {
// -----------------------------------------------------------------------------
Ana2Model::Ana2Model(const Geometry & resol, const eckit::Configuration & conf):
  geom_(new Geometry(resol)),
  rotate_(rotate(conf)),
  changegrid_(changegrid(conf))
{
  Log::trace() << "Ana2Model::Ana2Model start" << std::endl;
  Log::trace() << "Ana2Model::Ana2Model rotate:" << rotate_ << std::endl;
  Log::trace() << "Ana2Model::Ana2Model change grid:" << changegrid_ << std::endl;
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
  xm.rotate2grid();
  xa.agrid2uvgrid();
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
  xa.uvgrid2agrid();
  xa.rotate2north();
  xa.validTime() = xm.validTime();
  Log::trace() << "Ana2Model::changeVarInverse done" << xa << std::endl;
}
// -----------------------------------------------------------------------------
void Ana2Model::print(std::ostream & os) const {
  os << "Ana2Model";
}
// -----------------------------------------------------------------------------
}  // namespace soca
