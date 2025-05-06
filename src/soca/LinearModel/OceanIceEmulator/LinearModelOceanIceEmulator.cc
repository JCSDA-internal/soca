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
static oops::interface::LinearModelMaker<Traits, LinearModelOceanIceEmulator>
        makermodel_("LinearModelOceanIceEmulator");

  // -----------------------------------------------------------------------------
    LinearModelOceanIceEmulator::LinearModelOceanIceEmulator(const Geometry & resol,
                                        const eckit::Configuration & model)
        : tstep_(0), geom_(resol), traj_(), lrmodel_(geom_, model), vars_(model, "model variables")
    {
        Log::debug() << "------------ LinearModelOceanIceEmulator config: " << model << std::endl;
        tstep_ = util::Duration(model.getString("tstep"));
    }
    // -----------------------------------------------------------------------------
    LinearModelOceanIceEmulator::~LinearModelOceanIceEmulator() {
        Log::debug() << "------------ LinearModelOceanIceEmulator destructor" << std::endl;
        for (trajICst jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
            delete jtra->second;
        }
        Log::debug() << "------------ LinearModelOceanIceEmulator destructor done" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::initializeTL(Increment & dx) const {
        Log::debug() << "------------ LinearModelOceanIceEmulator::initializeTL" << std::endl;
        Log::debug() << "------------ dx:" << dx << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::initializeAD(Increment & dx) const {
        Log::debug() << "------------ LinearModelOceanIceEmulator::initializeAD" << std::endl;
        Log::debug() << "------------ dx:" << dx << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::stepTL(Increment & dx,
                                             const ModelBiasIncrement & bias) const {
        const ModelTrajectory * traj = this->getTrajectory(dx.validTime());
        oops::Log::debug() << "------------ Traj in TL" <<  traj << std::endl;
        dx.validTime() += tstep_;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::stepAD(Increment & dx, ModelBiasIncrement & bias) const {
        dx.validTime() -= tstep_;
        const ModelTrajectory * traj = this->getTrajectory(dx.validTime());
        oops::Log::debug() << "------------ Traj in AD" <<  traj << std::endl;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::setTrajectory(const State & xx,
                                                    State & xlr, const ModelBias & bias) {
        ASSERT(traj_.find(xx.validTime()) == traj_.end());
        ModelTrajectory * traj = new ModelTrajectory();
        traj->set(xlr);
        traj_[xx.validTime()] = traj;
    }
    // -----------------------------------------------------------------------------
    const ModelTrajectory * LinearModelOceanIceEmulator::getTrajectory(
                                                            const util::DateTime & tt) const {
      ASSERT(traj_.size() > 0);
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
