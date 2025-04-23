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
#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"

namespace soca {
// Define the Feed Forward Neural Net model
class OceanIceFFNN : public torch::nn::Module {
 public:
  OceanIceFFNN(const eckit::Configuration & config) {
    int inputSize = config.getInt("aimodel.inputSize");
    int hiddenSize = config.getInt("aimodel.hiddenSize");
    int outputSize = config.getInt("aimodel.outputSize");
    weightsFileName_ = config.getString("aimodel.load weights");
    normFileName_ = config.getString("aimodel.load normalization", "");

    oops::Log::info() << "Starting OceanIceFFNN constructor: "
                       << inputSize << outputSize << hiddenSize << std::endl;

    // Define the layers.
    fc1 = register_module("fc1", torch::nn::Linear(inputSize, hiddenSize));
    fc2 = register_module("fc2", torch::nn::Linear(hiddenSize, outputSize));

    // Register mean and std as buffers
    inputMean_ = register_buffer("input_mean", torch::full({inputSize}, 0.0));
    inputStd_ = register_buffer("input_std", torch::full({inputSize}, 1.0));
    oops::Log::trace() << "End OceanIceFFNN constructor" << std::endl;
  }

  // Initialize normalization
  void initNorm(torch::Tensor mean, torch::Tensor stdDev) {
    inputMean_ = mean;
    inputStd_ = stdDev;
  }

  void saveNorm(const std::string modelFileName) {
    std::vector<torch::Tensor> moments = {this->inputMean_, this->inputStd_};
    torch::save(moments, normFileName_);
  }

  void loadNorm() {
    std::vector<torch::Tensor> moments;
    torch::load(moments, normFileName_);
    this->inputMean_ = moments[0];
    this->inputStd_ = moments[1];
  }

  void initWeights() {
    // Initialize weights of the linear layers
    torch::nn::init::normal_(fc1->weight, 0.0, 0.5);  // Mean 0, StdDev 0.5
    torch::nn::init::normal_(fc2->weight, 0.0, 0.5);  // Mean 0, StdDev 0.5
  }

  void saveWeights() {
    std::vector<torch::Tensor> weights = {fc1->weight, fc2->weight};
    torch::save(weights, weightsFileName_);
  }

  void loadWeights() {
    // Load weights from file
    std::vector<torch::Tensor> weights;
    torch::load(weights, weightsFileName_);

    // Assign loaded weights to the layers
    fc1->weight = weights[0];
    fc2->weight = weights[1];
  }

  // Implement the forward pass with ReLU activation
  torch::Tensor forward(torch::Tensor x) {
    x = (x - inputMean_) / inputStd_;
    x = fc1(x);
    x = fc2(x);
    return x;
  }

  torch::Tensor forward_tlm(torch::Tensor x, torch::Tensor dx) {
    // Normalize the perturbation
    auto dx_norm = dx / inputStd_.detach();

    // TLM of the linear layers (no activation functions in your model)
    auto dx1 = fc1->weight.matmul(dx_norm.transpose(0, 1)).transpose(0, 1);  // Input → Hidden
    auto dx2 = fc2->weight.matmul(dx1.transpose(0, 1)).transpose(0, 1); // Hidden → Output
    return dx2;
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
  torch::Tensor inputMean_;
  torch::Tensor inputStd_;

  // Define the model name
  std::string weightsFileName_;
  std::string normFileName_;
};;
}  // namespace soca
