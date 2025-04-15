/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <filesystem>
#include <memory>
#include <utility>
# include <vector>
#include <string>

#include "torch/serialize.h"
#include "torch/torch.h"
#include "oops/util/Logger.h"

namespace soca {
// Define the Feed Forward Neural Net model
class OceanIceFFNN : public torch::nn::Module {
 public:
  OceanIceFFNN(const eckit::Configuration & config) {
    int inputSize = config.getInt("aimodel.inputSize");
    int hiddenSize = config.getInt("aimodel.hiddenSize");
    int outputSize = config.getInt("aimodel.outputSize");
    oops::Log::info() << "Starting OceanIceFFNN constructor: "
                       << inputSize << outputSize << hiddenSize << std::endl;

    // Define the layers.
    fc1 = register_module("fc1", torch::nn::Linear(inputSize, hiddenSize));
    fc2 = register_module("fc2", torch::nn::Linear(hiddenSize, outputSize));

    // Register mean and std as buffers
    inputMean = register_buffer("input_mean", torch::full({inputSize}, 0.0));
    inputStd = register_buffer("input_std", torch::full({inputSize}, 1.0));
    oops::Log::trace() << "End OceanIceFFNN constructor" << std::endl;
  }

  // Initialize normalization
  void initNorm(torch::Tensor mean, torch::Tensor stdDev) {
    inputMean = mean;
    inputStd = stdDev;
  }

  void saveNorm(const std::string modelFileName) {
    std::filesystem::path filePath(modelFileName);
    auto path = filePath.parent_path();
    auto fileName = filePath.filename();

    std::vector<torch::Tensor> moments = {this->inputMean, this->inputStd};
    torch::save(moments, path.string() + "/normalization." + fileName.string());
  }

  void loadNorm(const std::string modelFileName) {
    std::filesystem::path filePath(modelFileName);
    auto path = filePath.parent_path();
    auto fileName = filePath.filename();

    std::vector<torch::Tensor> moments;
    torch::load(moments, path.string() + "/normalization." + fileName.string());
    this->inputMean = moments[0];
    this->inputStd = moments[1];
  }

  void initWeights() {
    torch::nn::init::xavier_normal_(fc1->weight);
    torch::nn::init::xavier_normal_(fc2->weight);
  }

  // Implement the forward pass
  torch::Tensor forward(torch::Tensor x) {
    x = (x - inputMean) / inputStd;
    x = fc1(x);
    x = fc2(x);
    return x;
  }

torch::Tensor forward_tlm(torch::Tensor x, torch::Tensor dx) {
    // Ensure x requires gradients
    x = x.clone().detach().requires_grad_(true);

    // Compute the forward pass
    auto y = this->forward(x);

    // Compute the Jacobian-vector product (JVP)
    auto jvp = torch::autograd::grad(
        /* outputs */ torch::autograd::variable_list{y},  // Wrap y in a variable_list
        /* inputs */ torch::autograd::variable_list{x},  // Wrap x in a variable_list
        /* grad_outputs */ torch::autograd::variable_list{dx},  // Wrap dx in a variable_list
        /* retain_graph */ true,
        /* create_graph */ true
    );

    // Return the JVP
    return jvp[0];
}

torch::Tensor forward_ad(torch::Tensor x, torch::Tensor dy) {
    // Ensure x requires gradients
    x = x.clone().detach().requires_grad_(true);

    // Forward pass to get y
    auto y = this->forward(x);

    // Compute vector-Jacobian product (VJP), i.e., adjoint
    auto vjp = torch::autograd::grad(
        /* outputs */ torch::autograd::variable_list{y},
        /* inputs */ torch::autograd::variable_list{x},
        /* grad_outputs */ torch::autograd::variable_list{dy},
        /* retain_graph */ true,
        /* create_graph */ false  // No need for higher-order derivatives
    );

    // Return the VJP (i.e., gradient wrt input)
    return vjp[0];
}

 private:
  // Define the layers.
  torch::nn::Linear fc1{nullptr};
  torch::nn::Linear fc2{nullptr};
  torch::Tensor inputMean;
  torch::Tensor inputStd;
};
}  // namespace soca
