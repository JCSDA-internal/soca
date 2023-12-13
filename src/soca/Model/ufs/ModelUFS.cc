/*
 * (C) Copyright 2020 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ModelUFS.interface.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Model/ufs/ModelUFS.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"

namespace soca {
// -------------------------------------------------------------------------------------------------
// Options taken by ModelUFS
class ModelUFSParameters : public oops::ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(ModelUFSParameters, ModelParametersBase)

 public:
  oops::RequiredParameter<oops::Variables> modelVariables{ "model variables", this};
  oops::RequiredParameter<util::Duration> tstep{ "tstep", this};
  oops::RequiredParameter<std::string> ufsRunDirectory{ "ufs_run_directory", this};
};
// -------------------------------------------------------------------------------------------------
static oops::interface::ModelMaker<Traits, ModelUFS> makermodel_("UFS");
// -------------------------------------------------------------------------------------------------
ModelUFS::ModelUFS(const Geometry & resol, const eckit::Configuration & config)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_() {
  char tmpdir_[10000];
  oops::Log::trace() << "ModelUFS::ModelUFS starting" << std::endl;
  getcwd(tmpdir_, 10000);

  ModelUFSParameters params;
  params.deserialize(config);
  vars_ = params.modelVariables;
  tstep_ = params.tstep;
  strcpy(ufsdir_, params.ufsRunDirectory.value().c_str());
  chdir(ufsdir_);
  soca_ufs_create_f90(keyConfig_, config, geom_.toFortran());

  chdir(tmpdir_);
  oops::Log::trace() << "ModelUFS::ModelUFS done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelUFS::~ModelUFS() {
  oops::Log::trace() << "ModelUFS::~ModelUFS starting" << std::endl;
  soca_ufs_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelUFS::~ModelUFS done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::initialize(State & xx) const {
  oops::Log::trace() << "ModelUFS::initialize starting" << std::endl;
  oops::Log::trace() << "ModelUFS::cd to " << ufsdir_ << std::endl;

  chdir(ufsdir_);
  util::DateTime * dtp = &xx.validTime();
  soca_ufs_initialize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::trace() << "ModelUFS::initialize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::step(State & xx, const ModelBias &) const
{
  oops::Log::trace() << "ModelUFS::step starting" << std::endl;
  oops::Log::trace() << "ModelUFS::cd to " << ufsdir_ << std::endl;
  chdir(ufsdir_);

  util::DateTime start = xx.validTime();
  util::DateTime * dtp1 = &start;
  oops::Log::trace() << "Model start time is " << xx.validTime() << std::endl;
  oops::Log::trace() << "Forecast time step is " << tstep_ << std::endl;
  xx.validTime() += tstep_;
  util::DateTime * dtp2 = &xx.validTime();
  soca_ufs_step_f90(keyConfig_, xx.toFortran(), &dtp1, &dtp2);
  oops::Log::trace() << "ModelUFS::step done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::finalize(State & xx) const {
  oops::Log::trace() << "ModelUFS::finalize starting" << std::endl;
  util::DateTime * dtp = &xx.validTime();
  soca_ufs_finalize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::trace() << "ModelUFS::finalize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
int ModelUFS::saveTrajectory(State & xx, const ModelBias &) const {
//  ABORT("Model:UFS should not be used for the trajectory");
  return 0;
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::print(std::ostream & os) const {
  os << "ModelUFS::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace soca
