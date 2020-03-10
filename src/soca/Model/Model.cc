/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "soca/Traits.h"

#include "soca/Fields/Fields.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Model/Model.h"
#include "soca/Model/ModelFortran.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  static oops::ModelMaker<Traits, Model> makermodel_("SOCA");
  // -----------------------------------------------------------------------------
  Model::Model(const Geometry & resol, const eckit::Configuration & model)
    : tstep_(0), geom_(new Geometry(resol)), vars_(model),
      setup_mom6_(true)
  {
    Log::trace() << "Model::Model" << std::endl;
    Log::trace() << "Model vars: " << vars_ << std::endl;
    tstep_ = util::Duration(model.getString("tstep"));
    setup_mom6_ = model.getBool("setup_mom6", true);
    if (setup_mom6_)
      {
        soca_setup_f90(ftn_, &model, geom_->toFortran());
      }
        Log::trace() << "Model created" << std::endl;
  }
  // -----------------------------------------------------------------------------
  Model::~Model() {
    if (setup_mom6_)
      {
        soca_delete_f90(ftn_);
      }
    Log::trace() << "Model destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::initialize(State & xx) const {
    ASSERT(xx.fields().isForModel(true));
    soca_initialize_integration_f90(ftn_, xx.fields().toFortran());
    Log::debug() << "Model::initialize" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  void Model::step(State & xx, const ModelBias &) const {
    ASSERT(xx.fields().isForModel(true));
    Log::trace() << "Model::Time: " << xx.validTime() << std::endl;
    soca_propagate_f90(ftn_, xx.fields().toFortran(), &xx.validTime());
    xx.validTime() += tstep_;
  }
  // -----------------------------------------------------------------------------
  void Model::finalize(State & xx) const {
    ASSERT(xx.fields().isForModel(true));
    soca_finalize_integration_f90(ftn_, xx.fields().toFortran());
    Log::debug() << "Model::finalize" << xx.fields() << std::endl;
  }
  // -----------------------------------------------------------------------------
  int Model::saveTrajectory(State & xx, const ModelBias &) const {
    ASSERT(xx.fields().isForModel(true));
    int ftraj = 0;
    xx.validTime() += tstep_;
    return ftraj;
  }
  // -----------------------------------------------------------------------------
  void Model::print(std::ostream & os) const {
    os << "Model::print not implemented";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
