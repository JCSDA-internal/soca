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

// Define the Feed Forward Neural Net model
struct IceNet : torch::nn::Module {
  IceNet(int inputSize, int hiddenSize, int outputSize, int kernelSize = 1, int stride = 1) {
    oops::Log::trace() << "Starting IceNet constructor: "
                       << inputSize << outputSize << hiddenSize << std::endl;
    // Define Batch Normalization layer
    // batch_norm = register_module("batch_norm", torch::nn::BatchNorm1d(1, inputSize));

    // Define the layers.
    fc1 = register_module("fc1", torch::nn::Linear(inputSize, hiddenSize));
    fc2 = register_module("fc2", torch::nn::Linear(hiddenSize, outputSize));

    // Register mean and std as buffers
    inputMean = register_buffer("input_mean", torch::full({inputSize}, 0.0));
    inputStd = register_buffer("input_std", torch::full({inputSize}, 1.0));
    oops::Log::trace() << "End IceNet constructor" << std::endl;
  }

  // Initialize normalization
  void initNorm(torch::Tensor mean, torch::Tensor std) {
    inputMean = mean;
    inputStd = std;
  }

  void saveNorm(const std::string modelFileName) {
    // construct normalization file name
    // TODO(G): Move as a data member
    std::filesystem::path filePath(modelFileName);
    auto path = filePath.parent_path();
    auto fileName = filePath.filename();

    // Save 1st and 2nd moments
    std::vector<torch::Tensor> moments = {this->inputMean, this->inputStd};
    torch::save(moments, path.string() + "/normalization." + fileName.string());
  }

  void loadNorm(const std::string modelFileName) {
    // construct normalization file name
    std::filesystem::path filePath(modelFileName);
    auto path = filePath.parent_path();
    auto fileName = filePath.filename();

    // Load 1st and 2nd moments
    std::vector<torch::Tensor> moments;
    torch::load(moments, path.string() + "/normalization." + fileName.string());
    this->inputMean = moments[0];
    this->inputStd = moments[1];
  }

  void initWeights() {
    // Xavier initialization for the first two layers
    torch::nn::init::xavier_normal_(fc1->weight);
    torch::nn::init::xavier_normal_(fc2->weight);
  }

  // void load(const std::string modelFileName) {
  // torch::load(*this, modelFileName);
    /*
    for (const auto& pair : this->named_parameters()) {
      std::cout << "Parameter name: " << pair.key()
                << ", Size: " << pair.value().sizes() << std::endl;
    }
    for (const auto& pair : this->named_buffers()) {
      std::cout << "Buffer name: " << pair.key()
                << ", Size: " << pair.value().sizes() << std::endl;
      std::cout << "       values: " << pair.value() << std::endl;
      }*/
  //}

  // Implement the forward pass
  torch::Tensor forward(torch::Tensor x) {
    // Normalize the input
    x = (x - inputMean) / inputStd;
    x = fc1(x);
    x = torch::sigmoid(fc2(x));
    return x;
  }

  // Compute the Jacobian (dout/dx)
  torch::Tensor jac(torch::Tensor x) {
    auto xp = torch::from_blob(x.data_ptr(), {x.size(0)}, torch::requires_grad());
    auto y = this->forward(xp);
    y.backward();
    return xp.grad();
  }

  torch::Tensor jacNorm(torch::Tensor input) {
    torch::Tensor frobeniusNorm = torch::tensor(0.0);
    return frobeniusNorm;
  }
  // Define the layers.
  torch::nn::Linear fc1{nullptr};
  torch::nn::Linear fc2{nullptr};
  // torch::nn::BatchNorm1d batch_norm{nullptr};
  torch::Tensor inputMean;
  torch::Tensor inputStd;
};
