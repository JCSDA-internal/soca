/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Model/Model.h"

#include "oops/util/Logger.h"
#include "src/ModelBias.h"
#include "src/Fields/Fields.h"
#include "src/Fortran.h"
#include "src/Geometry/Geometry.h"
#include "src/State/State.h"
#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  Model::Model(const Geometry & resol, const eckit::Configuration & model)
    : keyConfig_(0), tstep_(0), geom_(resol)
  {
    Log::trace() << "Model::Model" << std::endl;
    tstep_ = util::Duration(model.getString("tstep"));
    const eckit::Configuration * configc = &model;
    soca_setup_f90(&configc, geom_.toFortran(), keyConfig_);
    Log::trace() << "Model created" << std::endl;
  }
  // -----------------------------------------------------------------------------
  Model::~Model() {
    soca_delete_f90(keyConfig_);
    Log::trace() << "Model destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::initialize(State & xx) const {
    xx.activateModel();
    ASSERT(xx.fields().isForModel(true));
    soca_prepare_integration_f90(keyConfig_, xx.fields().toFortran());
    Log::debug() << "Model::initialize" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::step(State & xx, const ModelBias &) const {
    ASSERT(xx.fields().isForModel(true));
    Log::debug() << "Model::step fields in" << xx.fields() << std::endl;
    soca_propagate_f90(keyConfig_, xx.fields().toFortran());
    xx.validTime() += tstep_;
    Log::debug() << "Model::step fields out" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::finalize(State & xx) const {
    xx.deactivateModel();
    Log::debug() << "Model::finalize" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  int Model::saveTrajectory(State & xx, const ModelBias &) const {
    ASSERT(xx.fields().isForModel(true));
    int ftraj = 0;
    Log::debug() << "Model::saveTrajectory fields in" << xx.fields() << std::endl;
    //soca_prop_traj_f90(keyConfig_, xx.fields().toFortran(), ftraj);
    xx.validTime() += tstep_;
    ASSERT(ftraj != 0);
    Log::debug() << "Model::saveTrajectory fields out" << xx.fields() << std::endl;
    return ftraj;
  }
  // -----------------------------------------------------------------------------
  void Model::print(std::ostream & os) const {
    os << "Model::print not implemented";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
