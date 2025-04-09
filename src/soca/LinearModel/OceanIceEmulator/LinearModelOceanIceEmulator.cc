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
#include "soca/LinearModel/OceanIceEmulator/LinearModelOceanIceEmulator.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"
#include "soca/Increment/Increment.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace soca {
  // -----------------------------------------------------------------------------
  static oops::interface::LinearModelMaker<Traits, LinearModelOceanIceEmulator> makermodel_("LinearModelOceanIceEmulator");

  // -----------------------------------------------------------------------------
    LinearModelOceanIceEmulator::LinearModelOceanIceEmulator(const Geometry & resol,
                                        const eckit::Configuration & model)
        : tstep_(0), geom_(resol), traj_(), lrmodel_(geom_, model), vars_(model, "model variables")
    {
        Log::trace() << "------------ LinearModelOceanIceEmulator::LinearModelOceanIceEmulator" << std::endl;
        Log::trace() << "------------ LinearModelOceanIceEmulator config: " << model << std::endl;
        tstep_ = util::Duration(model.getString("tstep"));
    }
    // -----------------------------------------------------------------------------
    LinearModelOceanIceEmulator::~LinearModelOceanIceEmulator() {
        Log::trace() << "------------ LinearModelOceanIceEmulator destructed" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::initializeTL(Increment & dx) const {
        Log::debug() << "------------ LinearModelOceanIceEmulator::initializeTL" << std::endl;
        Log::debug() << "------------ dx:" << dx << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::initializeAD(Increment & dx) const {
        Log::debug() << "------------ LinearModelOceanIceEmulator::initializeAD" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::stepTL(Increment & dx, const ModelBiasIncrement & bias) const {
        Log::trace() << "------------ LinearModelOceanIceEmulator::Time: " << dx.validTime() << std::endl;
        //util::DateTime * modeldate = &dx.validTime();
        const ModelTrajectory * traj = this->getTrajectory(dx.validTime());
        dx.validTime() += tstep_;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::stepAD(Increment & dx, ModelBiasIncrement & bias) const {
        Log::trace() << "------------ LinearModelOceanIceEmulator::Time: " << dx.validTime() << std::endl;
        //util::DateTime * modeldate = &dx.validTime();
        dx.validTime() -= tstep_;
        const ModelTrajectory * traj = this->getTrajectory(dx.validTime());
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::setTrajectory(const State & xx, State & xlr, const ModelBias & bias) {
        Log::trace() << "------------ LinearModelOceanIceEmulator::setTrajectory" << std::endl;
        Log::trace() << "xlr: " << xlr << std::endl;
        ASSERT(traj_.find(xx.validTime()) == traj_.end());
        ModelTrajectory * traj = new ModelTrajectory();
        Log::info() << "$$$$$$$$$$$$$$$$$$$$$$$$$ set traj time: " << xx.validTime() << std::endl;
        //lrmodel_.step(xlr, bias);
        traj->set(xlr);
        Log::info() << "------------ LinearModelOceanIceEmulator::setTrajectory " << xx << std::endl;
        traj_[xx.validTime()] = traj;
    }
    // -----------------------------------------------------------------------------
    const ModelTrajectory * LinearModelOceanIceEmulator::getTrajectory(const util::DateTime & tt) const {
      ASSERT(traj_.begin()->first <= tt);
      ASSERT(traj_.rbegin()->first >= tt);
      trajICst itra = traj_.lower_bound(tt);
      return itra->second;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::finalizeTL(Increment & dx) const {
        Log::debug() << "LinearModelOceanIceEmulator::finalizeTL" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::finalizeAD(Increment & dx) const {
        Log::debug() << "LinearModelOceanIceEmulator::finalizeAD" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::print(std::ostream & os) const {
        Log::info() << "LinearModelOceanIceEmulator::print not implemented" << std::endl;
    }
    // -----------------------------------------------------------------------------
}  // namespace soca
