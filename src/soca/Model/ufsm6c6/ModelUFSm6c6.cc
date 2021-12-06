/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "soca/Traits.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Model/ufsm6c6/ModelUFSm6c6.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  static oops::interface::ModelMaker<soca::Traits, soca::ModelUFSm6c6>
               makermodel_("UFSm6c6");
  // -----------------------------------------------------------------------------
  ModelUFSm6c6::ModelUFSm6c6(const Geometry & resol,
                             const eckit::Configuration & model)
    : keyConfig_(0),
      tstep_(0),
      geom_(new Geometry(resol)),
      vars_(model, "model variables")
  {
    Log::trace() << "ModelUFSm6c6::ModelUFSm6c6" << std::endl;
    Log::trace() << "ModelUFSm6c6 vars: " << vars_ << std::endl;
    tstep_ = util::Duration(model.getString("tstep"));
  }
  // -----------------------------------------------------------------------------
  ModelUFSm6c6::~ModelUFSm6c6() {
    Log::trace() << "ModelUFSm6c6 destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void ModelUFSm6c6::initialize(State & xx) const {
    Log::debug() << "ModelUFSm6c6::initialize" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void ModelUFSm6c6::step(State & xx, const ModelBias &) const {
    Log::trace() << "ModelUFSm6c6::Time: " << xx.validTime() << std::endl;
    util::DateTime * modeldate = &xx.validTime();
    xx.validTime() += tstep_;
  }
  // -----------------------------------------------------------------------------
  void ModelUFSm6c6::finalize(State & xx) const {
    Log::debug() << "ModelUFSm6c6::finalize" << std::endl;
  }
  // -----------------------------------------------------------------------------
  int ModelUFSm6c6::saveTrajectory(State & xx, const ModelBias &) const {
    int ftraj = 0;
    xx.validTime() += tstep_;
    return ftraj;
  }
  // -----------------------------------------------------------------------------
  void ModelUFSm6c6::print(std::ostream & os) const {
    os << "ModelUFSm6c6::print not implemented";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
