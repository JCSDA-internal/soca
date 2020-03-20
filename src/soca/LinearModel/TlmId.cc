/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <ostream>
#include <vector>

#include "soca/Traits.h"

#include "soca/Fortran.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/LinearModel/TlmId.h"
#include "soca/ModelBias/ModelBiasIncrement.h"
#include "soca/State/State.h"

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
// -----------------------------------------------------------------------------
static oops::LinearModelMaker<Traits, TlmId> makerIdTLM_("IdTLM");
// -----------------------------------------------------------------------------
TlmId::TlmId(const Geometry & resol, const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol), linvars_(tlConf)
{
  tstep_ = util::Duration(tlConf.getString("tstep"));
  Log::trace() << "TlmId created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmId::~TlmId() {
  Log::trace() << "TlmId destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::setTrajectory(const State &, State &, const ModelBias &) {}
// -----------------------------------------------------------------------------
void TlmId::initializeTL(Increment & dx) const {}
// -----------------------------------------------------------------------------
void TlmId::stepTL(Increment & dx, const ModelBiasIncrement &) const {
  dx.updateTime(tstep_);
}
// -----------------------------------------------------------------------------
void TlmId::finalizeTL(Increment & dx) const {}
// -----------------------------------------------------------------------------
void TlmId::initializeAD(Increment & dx) const {}
// -----------------------------------------------------------------------------
void TlmId::stepAD(Increment & dx, ModelBiasIncrement &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmId::finalizeAD(Increment & dx) const {
  Log::debug() << "TlmId::finalizeAD" << dx << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::print(std::ostream & os) const {
  os << " IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace soca
