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
#include "soca/Model/OceanIceEmulator/OceanIceFFNN.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/base/GeometryData.h"

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
      geomData_(resol.functionSpace(), resol.fields(),
                resol.levelsAreTopDown(), resol.getComm()),
      vars_(model, "model variables"),
      aimodel_(model)
  {
    // Get the time step from the model configuration
    tstep_ = util::Duration(model.getString("tstep"));

    // Load the AI model weights
    aimodel_.loadWeights();
    Log::debug() << "------------ Loaded AI model weights in NL model: " << std::endl;
    for (const auto& tensor : aimodel_.parameters()) {
      Log::debug() << tensor << std::endl;
    }
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
    Log::debug() << "------------ ModelOceanIceEmulator::Time: " << xx.validTime() << std::endl;

    // Create a vector of field views
    std::vector<atlas::array::ArrayView<double, 2>> x_v;
    std::vector<std::string> x_names;

    // Populate the vector with views to each field
    for (auto & field : xx.fieldSet()) {
      oops::Log::debug() << "------------ field: " << field.name() << std::endl;
      x_v.push_back(atlas::array::make_view<double, 2>(field));
      x_names.push_back(field.name());
    }

    // mask and ghost cells
    auto fs = geomData_.functionSpace();
    const auto & ghostView = atlas::array::make_view<int, 1>(fs.ghost());
    const auto & maskView = atlas::array::make_view<double, 2>(geomData_.getField("mask_h"));

    // Print the weights of the aimodel
    for (const auto& tensor : aimodel_.parameters()) {
        Log::debug() << tensor << " ";
    }
    Log::debug() << std::endl;

    // Loop through nodes and apply the aimodel
    auto nodes = x_v[0].shape(0);
    for (size_t jnode = 0; jnode < nodes; ++jnode) {
      // Skip land points and ghost points
      if (maskView(jnode, 0) == 0 || ghostView(jnode)) continue;  // skip land/ghost points
      //if (maskView(jnode, 0) == 0 || ghostView(jnode)) {
      //    for (size_t j = 0; j < x_v.size(); ++j) {
      //      x_v[j](jnode, 0) = 0.0;
      //  }
      //  continue;  // skip land/ghost points
      //}

      // Prepare the input data
      std::vector<double> inputData;
      for (size_t j = 0; j < x_v.size(); ++j) {
        inputData.push_back(x_v[j](jnode, 0));
      }

      // Convert inputData to a tensor
      torch::Tensor inputTensor = torch::tensor(inputData).reshape({1, static_cast<int64_t>(inputData.size())});
      auto torch_x_out = aimodel_.forward(inputTensor);

      // Update the state
      for (size_t j = 0; j < x_v.size(); ++j) {
        x_v[j](jnode, 0) = torch_x_out[0][j].item<double>();
      }
    }
    // Update halo points
    //for (auto & field : xx.fieldSet()) {
    //  field.haloExchange();
    //}
    // Update the valid time
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
