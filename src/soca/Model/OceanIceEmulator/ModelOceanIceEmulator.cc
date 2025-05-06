/*
* (C) Copyright 2025 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <vector>

#include "soca/Traits.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Model/OceanIceEmulator/ModelOceanIceEmulator.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  static oops::interface::ModelMaker<soca::Traits, soca::ModelOceanIceEmulator>
      makermodel_("ModelOceanIceEmulator");

  // -----------------------------------------------------------------------------
  ModelOceanIceEmulator::ModelOceanIceEmulator(const Geometry & resol,
                             const eckit::Configuration & model)
    : tstep_(0),
      geom_(resol),
      vars_(model, "model variables")
  {
    Log::trace() << "------------ ModelOceanIceEmulator::ModelOceanIceEmulator" << std::endl;
    Log::trace() << "------------ ModelOceanIceEmulator vars: " << vars_ << std::endl;
    tstep_ = util::Duration(model.getString("tstep"));
  }
  // -----------------------------------------------------------------------------
  ModelOceanIceEmulator::~ModelOceanIceEmulator() {
    Log::trace() << "------------ ModelOceanIceEmulator destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void ModelOceanIceEmulator::initialize(State & xx) const {
    Log::debug() << "------------ ModelOceanIceEmulator::initialize" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void ModelOceanIceEmulator::step(State & xx, const ModelBias &) const {
    Log::trace() << "------------ ModelOceanIceEmulator::Time: " << xx.validTime() << std::endl;
    xx.validTime() += tstep_;
  }
  // -----------------------------------------------------------------------------
  void ModelOceanIceEmulator::finalize(State & xx) const {
    Log::debug() << "------------ ModelOceanIceEmulator::finalize" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void ModelOceanIceEmulator::print(std::ostream & os) const {
    os << "ModelOceanIceEmulator::print not implemented";
  }
  // -----------------------------------------------------------------------------
}  // namespace soca
