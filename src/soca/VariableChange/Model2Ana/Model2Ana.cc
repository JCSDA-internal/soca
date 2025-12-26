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
#include "oops/util/Timer.h"

#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"
#include "soca/VariableChange/Model2Ana/Model2Ana.h"

using oops::Log;

namespace soca {

// -----------------------------------------------------------------------------

static VariableChangeMaker<Model2Ana> makerVarChaA2M_("Model2Ana");

// -----------------------------------------------------------------------------
Model2Ana::Model2Ana(const Geometry & resol, const eckit::Configuration & conf)
  : uvars_(initRotate(conf, "u")), vvars_(initRotate(conf, "v")),
    interp_(initInterp(conf)), logvars_(initTrans(conf, "var"))
{
  Log::trace() << "Model2Ana::Model2Ana start" << std::endl;
  util::Timer timer("soca::Model2Ana", "Model2Ana");

  ASSERT(uvars_.size() == vvars_.size());
  Log::trace() << "Model2Ana::Model2Ana Rotating:"
               << " u = " << uvars_ << " v = " << vvars_ << std::endl;
  Log::trace() << "Model2Ana::Model2Ana Interpolate:"
               << interp_ << std::endl;
  Log::trace() << "Model2Ana::Log Transforming:"
               << " var = " << logvars_ << std::endl;
  Log::trace() << "Model2Ana::Model2Ana done" << std::endl;
}
// -----------------------------------------------------------------------------
Model2Ana::~Model2Ana() {
  oops::Log::trace() << "ChangeSOCA destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void Model2Ana::changeVar(const State & xm,
                          State & xa) const {
  Log::trace() << "Model2Ana::changeVar starting" <<xm <<
                        std::endl;
  util::Timer timer("soca::Model2Ana", "changeVar");

  util::DateTime * vtime = &xa.validTime();
  xa = xm;

  // Rotate from the logical grid to meridional/zonal
  xa.rotate2north(uvars_, vvars_);
  Log::trace() << "Model2Ana::changeVar rotated to meridional/zonal" << xa << std::endl;

  // If the interp switch is true, interpolate vectors to the h-grid
  if (interp_) {
    xa.tohgrid(uvars_, vvars_);
    Log::trace() << "Model2Ana::changeVar interpolated to h-grid" << xa << std::endl;
  }

  // Apply Log transform on the variables in logvars_
  xa.logtrans(logvars_);
  Log::trace() << "Model2Ana::changeVar applied Log transform" << xa << std::endl;

  xa.validTime() = xm.validTime();
  Log::trace() << "Model2Ana::changeVar done" << xa << std::endl;
}
// -----------------------------------------------------------------------------
void Model2Ana::changeVarInverse(const State & xa,
                                State & xm) const {
  oops::Log::trace() << "Model2Ana::changeVarInverse starting" << xa <<
                        std::endl;
  util::Timer timer("soca::Model2Ana", "changeVarInverse");

  util::DateTime * vtime = &xm.validTime();
  xm = xa;

  // Rotate from meridional/zonal to the logical grid
  xm.rotate2grid(uvars_, vvars_);
  Log::trace() << "Model2Ana::changeVarInverse rotated to logical grid" << xm << std::endl;

  // If the interp switch is true, interpolate vectors back to the c-grid
  if (interp_) {
    xm.tocgrid(uvars_, vvars_);
    Log::trace() << "Model2Ana::changeVarInverse interpolated vectors to c-grid" << xm << std::endl;
  }

  // exp transform the variables in logvars_
  xm.expontrans(logvars_);
  Log::trace() << "Model2Ana::changeVarInverse applied Log transform" << xm << std::endl;

  xm.validTime() = xa.validTime();
  Log::trace() << "Model2Ana::changeVarInverse done" << xm << std::endl;
}
// -----------------------------------------------------------------------------
void Model2Ana::print(std::ostream & os) const {
  os << "Model2Ana";
}
// -----------------------------------------------------------------------------
}  // namespace soca
