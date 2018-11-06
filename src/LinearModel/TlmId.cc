/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "src/LinearModel/TlmId.h"

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/Logger.h"
#include "src/ModelBiasIncrement.h"
#include "src/Fortran.h"
#include "src/Geometry/Geometry.h"
#include "src/Increment/Increment.h"
#include "src/State/State.h"
#include "src/Traits.h"
#include "oops/util/DateTime.h"
#include "oops/util/abor1_cpp.h"


using oops::Log;

namespace soca {
// -----------------------------------------------------------------------------
static oops::LinearModelMaker<Traits, TlmId> makerIdTLM_("IdTLM");
// -----------------------------------------------------------------------------
TlmId::TlmId(const Geometry & resol, const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol),
    linvars_(std::vector<std::string>{
          "cicen",
          "hicen",
          "hsnon",
          "tsfcn",
          "qsnon",
          "sicnk",
          "qicnk",
          "socn",
          "tocn",
          "ssh",
          "hocn"})
{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;
  soca_setup_f90(&configc, resol_.toFortran(), keyConfig_);

  Log::trace() << "TlmId created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmId::~TlmId() {
  soca_delete_f90(keyConfig_);
  Log::trace() << "TlmId destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::setTrajectory(const State &, State &, const ModelBias &) {}
// -----------------------------------------------------------------------------
void TlmId::initializeTL(Increment & dx) const {
  dx.activateModel();
  ASSERT(dx.fields().isForModel(false));
  // soca_prepare_integration_tl_f90(keyConfig_, dx.fields().toFortran());
  Log::debug() << "TlmId::initializeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::stepTL(Increment & dx, const ModelBiasIncrement &) const {
  dx.updateTime(tstep_);
}
// -----------------------------------------------------------------------------
void TlmId::finalizeTL(Increment & dx) const {
  dx.deactivateModel();
  Log::debug() << "TlmId::finalizeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::initializeAD(Increment & dx) const {
  dx.activateModel();
  ASSERT(dx.fields().isForModel(false));
  Log::debug() << "TlmId::initializeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::stepAD(Increment & dx, ModelBiasIncrement &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmId::finalizeAD(Increment & dx) const {
  // soca_prepare_integration_ad_f90(keyConfig_, dx.fields().toFortran());
  dx.deactivateModel();
  Log::debug() << "TlmId::finalizeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::print(std::ostream & os) const {
  os << " IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace soca
