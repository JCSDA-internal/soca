#pragma once

#include <iostream>
#include <torch/torch.h>

#include "soca/Model/OceanIceEmulator/OceanIceFFNN.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/PostProcessor.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"


using soca::OceanIceFFNN;

namespace soca {
  double dot_local(torch::Tensor a, torch::Tensor b) {
    return torch::dot(a.flatten(), b.flatten()).item<double>();
  }

  void printTestResult(const std::string &testName, double err, double tol) {
      bool passed = err < tol;
      oops::Log::info() << "=== Test: " << testName << "\n";
      oops::Log::info() << "Error: " << err << "\n";
      oops::Log::info() << "Tolerance: " << tol << "\n";
      if (passed) {
          oops::Log::info() << "\033[1m\033[32m[PASS]\033[0m " << testName << "\n";
      } else {
          oops::Log::info() << "\033[1m\033[31m[FAIL]\033[0m " << testName << "\n";
      }
  }

  class test_OceanIceFFNN : public oops::Application {
    public:
    explicit test_OceanIceFFNN(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}
    static const std::string classname() {return "soca::test_OceanIceFFNN";}

    int execute(const eckit::Configuration & fullConfig) const {
      // Setup OceanIceEmulator
      oops::Log::debug() << "--- OceanIceEmulator ctor" << std::endl;
      OceanIceFFNN model(fullConfig);

      // Create model
      oops::Log::debug() << "--- OceanIceEmulator initWeights" << std::endl;
      model.initWeights();
      oops::Log::debug() << "--- OceanIceEmulator saveWeights" << std::endl;
      model.saveWeights();
      oops::Log::debug() << "--- OceanIceEmulator initNorm" << std::endl;
      model.initNorm(torch::zeros({2}), torch::ones({2}));

      // Generate random inputs
      oops::Log::debug() << "--- OceanIceEmulator gen random input" << std::endl;
      auto x = torch::randn({1, 2}, torch::dtype(torch::kDouble));
      auto dx = torch::randn_like(x, torch::dtype(torch::kDouble));
      auto dy = torch::randn({1, 2}, torch::dtype(torch::kDouble));

      // Forward base output
      oops::Log::debug() << "--- OceanIceEmulator forward" << std::endl;
      auto y = model.forward(x);

      // ---------- TLM Check: finite difference ----------
      oops::Log::debug() << "--- OceanIceEmulator TLM check: finite diff" << std::endl;
      double eps = 1e-5;
      auto y_fd = model.forward(x + eps*dx);
      auto dy_fd = (y_fd - y);
      auto dy_tlm = model.forward_tlm(x, eps*dx);
      auto diff_tlm = torch::max(torch::abs(dy_fd - dy_tlm)).item<double>();
      printTestResult("TLM finite difference test", diff_tlm, 1e-6);

      // Check linearity of TLM
      oops::Log::debug() << "--- OceanIceEmulator TLM check: linearity" << std::endl;
      double scale = 2.0;
      auto y_tlm1 = model.forward_tlm(x, scale*dx);
      auto y_tlm2 = scale * model.forward_tlm(x, dx);
      auto diff_tlm2 = torch::max(torch::abs(y_tlm1 - y_tlm2)).item<double>();
      printTestResult("TLM linearity test", diff_tlm2, 1e-6);

      // ---------- AD Check: dot product test ----------
      auto tlm = model.forward_tlm(x, dx);
      auto adj = model.forward_ad(x, dy);

      double dot1 = dot_local(tlm, dy);
      double dot2 = dot_local(dx, adj);
      double relative_diff = std::abs(dot1 - dot2) / std::abs(dot1);
      printTestResult("Adjoint dot product test", relative_diff, 1e-6);

      return 0;
    }
   private:
    std::string appname() const {
      return "soca::test_OceanIceFFNN<";
    }
  };
}
