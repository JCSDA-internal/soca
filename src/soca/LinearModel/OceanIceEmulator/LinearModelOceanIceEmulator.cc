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
        : tstep_(0),
          geom_(resol),
          geomData_(resol.functionSpace(), resol.fields(),
                    resol.levelsAreTopDown(), resol.getComm()),
          ghostView_(atlas::array::make_view<int, 1>(geomData_.functionSpace().ghost())),
          maskView_(atlas::array::make_view<double, 2>(geomData_.getField("mask_h"))),
          traj_(),
          aimodel_(model),
          vars_(model, "model variables")
    {
        // Get the time step from the model configuration
        tstep_ = util::Duration(model.getString("tstep"));

        // Load the AI model weights
        aimodel_.loadWeights();

        // Ensure weights are detached
        for (auto &tensor : aimodel_.parameters()) {
            tensor = tensor.detach();  // Detach the tensor to avoid in-place operations
        }
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
    void LinearModelOceanIceEmulator::initializeTL(Increment & dx) const {}
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::initializeAD(Increment & dx) const {}
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::stepTL(Increment & dx,
                                             const ModelBiasIncrement & bias) const {
        // Get the trajectory
        auto [traj_v, traj_names] = this->getTrajectory(dx.validTime());

        // Create vectors of increment field views and names
        std::vector<atlas::array::ArrayView<double, 2>> dx_v;
        std::vector<std::string> dx_names;
        for (auto & field : dx.fieldSet()) {
          dx_v.push_back(atlas::array::make_view<double, 2>(field));
          dx_names.push_back(field.name());
        }

        // Loop through nodes and apply the linearized aimodel
        auto nodes = dx_v[0].shape(0);
        for (size_t jnode = 0; jnode < nodes; ++jnode) {
          // Skip land points and ghost points
          if (maskView_(jnode, 0) == 0 || ghostView_(jnode)) continue;
          //    for (size_t j = 0; j < dx_v.size(); ++j) {
          //      dx_v[j](jnode, 0) = 0.0;
          //  }
          //  continue;  // skip land/ghost points
          //}

          // Prepare the incrememnt input
          std::vector<double> inputData;
          for (size_t j = 0; j < dx_v.size(); ++j) {
            inputData.push_back(dx_v[j](jnode, 0));
          }

          // Prepare the trajectory input
          std::vector<double> inputDataTraj;
          for (size_t j = 0; j < traj_v.size(); ++j) {
            inputDataTraj.push_back(traj_v[j](jnode, 0));
          }

          // Convert inputData to a tensor
          torch::Tensor torch_dx = torch::tensor(inputData, torch::dtype(torch::kDouble))
                                       .reshape({1, static_cast<int64_t>(inputData.size())});
          torch::Tensor torch_traj = torch::tensor(inputDataTraj, torch::dtype(torch::kDouble))
                                         .reshape({1, static_cast<int64_t>(inputDataTraj.size())});
          auto torch_dx_out = aimodel_.forward_tlm(torch_traj, torch_dx);

          // Update the increment
          for (size_t j = 0; j < dx_v.size(); ++j) {
              dx_v[j](jnode, 0) = torch_dx_out[0][j].item<double>();
          }
        }
        // Update halo points
        //for (auto & field : dx.fieldSet()) {
        //    field.haloExchange();
        //}
        // Update the valid time
        dx.validTime() += tstep_;
    }
    // -----------------------------------------------------------------------------
    void LinearModelOceanIceEmulator::stepAD(Increment & dx, ModelBiasIncrement & bias) const {
        dx.validTime() -= tstep_;
        // Update halo points
        //for (auto & field : dx.fieldSet()) {
        //    field.haloExchange();
        //}

        // Get the trajectory
        auto [traj_v, traj_names] = this->getTrajectory(dx.validTime());

        // Create vectors of increment field views and names
        std::vector<atlas::array::ArrayView<double, 2>> dx_v;
        std::vector<std::string> dx_names;
        for (auto & field : dx.fieldSet()) {
          dx_v.push_back(atlas::array::make_view<double, 2>(field));
          dx_names.push_back(field.name());
        }
        // Loop through nodes and apply the linearized aimodel
        auto nodes = dx_v[0].shape(0);
        for (size_t jnode = 0; jnode < nodes; ++jnode) {
        //for (int jnode = nodes - 1; jnode >= 0; --jnode) {
          if (maskView_(jnode, 0) == 0 || ghostView_(jnode)) continue;  // skip land/ghost points
          //if (maskView_(jnode, 0) == 0 || ghostView_(jnode)) {
          //    for (size_t j = 0; j < dx_v.size(); ++j) {
          //      dx_v[j](jnode, 0) = 0.0;
          //  }
          //  continue;  // skip land/ghost points
          //}

          // Prepare the incrememnt input
          std::vector<double> inputData;
          for (size_t j = 0; j < dx_v.size(); ++j) {
            inputData.push_back(dx_v[j](jnode, 0));
          }

          // Prepare the trajectory input
          std::vector<double> inputDataTraj;
          for (size_t j = 0; j < traj_v.size(); ++j) {
            inputDataTraj.push_back(traj_v[j](jnode, 0));
          }

          // Convert inputData to a tensor
          torch::Tensor torch_dx = torch::tensor(inputData, torch::dtype(torch::kDouble))
                                       .reshape({1, static_cast<int64_t>(inputData.size())});
          torch::Tensor torch_traj = torch::tensor(inputDataTraj, torch::dtype(torch::kDouble))
                                         .reshape({1, static_cast<int64_t>(inputDataTraj.size())});
          auto torch_dx_out = aimodel_.forward_ad(torch_traj, torch_dx);

          // Update the increment
          for (size_t j = 0; j < dx_v.size(); ++j) {
              dx_v[j](jnode, 0) = torch_dx_out[0][j].item<double>();
          }
        }
        // Update halo points
        //for (auto & field : dx.fieldSet()) {
        //    field.haloExchange();
        //}
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

    std::pair<std::vector<atlas::array::ArrayView<double, 2>>,
              std::vector<std::string>>
    LinearModelOceanIceEmulator::getTrajectory(const util::DateTime & tt) const {
        ASSERT(traj_.size() > 0);  // Ensure the map is not empty

        // Find the first key greater than tt
        auto itra = traj_.upper_bound(tt);

        // If tt is earlier than the first key, throw an error
        ASSERT(itra != traj_.begin());

        // Move back to the trajectory valid for tt
        --itra;

        // Ensure that the trajectory is valid for tt
        ASSERT(itra->first <= tt);

        // Retrieve the ModelTrajectory
        const ModelTrajectory * traj = itra->second;

        // Retrieve the State valid for tt
        soca::State trajState = traj->getStateByDateTime(tt);
        //return traj->getStateByDateTime(tt);

        // Create trajectory views
        std::vector<atlas::array::ArrayView<double, 2>> traj_v;
        std::vector<std::string> traj_names;
        for (auto & field : trajState.fieldSet()) {
          traj_v.push_back(atlas::array::make_view<double, 2>(field));
          traj_names.push_back(field.name());
        }
        return std::make_pair(traj_v, traj_names);
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
