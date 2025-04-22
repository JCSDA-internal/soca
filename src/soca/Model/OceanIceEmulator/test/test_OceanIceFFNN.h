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

  void printTestResult(const std::string &testName, bool passed) {
      if (passed) {
          std::cout << "\033[1m\033[32m[PASS]\033[0m " << testName << "\n";
      } else {
          std::cout << "\033[1m\033[31m[FAIL]\033[0m " << testName << "\n";
      }
  }

  class test_OceanIceFFNN : public oops::Application {
    public:
    explicit test_OceanIceFFNN(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}
    static const std::string classname() {return "soca::test_OceanIceFFNN";}

    int execute(const eckit::Configuration & fullConfig) const {
      // Setup OceanIceEmulator
      OceanIceFFNN model(fullConfig);
      oops::Log::info() << "------------ OceanIceFFNN config: " << fullConfig << std::endl;

      // Create model
      model.initWeights();
      model.initNorm(torch::zeros({2}), torch::ones({2}));
      oops::Log::info() << "------------ OceanIceFFNN initialized" << std::endl;

      // Generate random inputs
      auto x = torch::randn({1, 2});
      auto dx = torch::randn_like(x);
      auto dy = torch::randn({1, 2});
      oops::Log::info() << "------------ OceanIceFFNN inputs: " << x << std::endl;
      oops::Log::info() << "------------ OceanIceFFNN dx: " << dx << std::endl;
      oops::Log::info() << "------------ OceanIceFFNN dy: " << dy << std::endl;

      // Forward base output
      auto y = model.forward(x);
      oops::Log::info() << "------------ OceanIceFFNN forward output: " << y << std::endl;

      // ---------- TLM Check: finite difference ----------
      double eps = 1e-5;
      auto y_fd = model.forward(x + eps*dx);
      auto dy_fd = (y_fd - y);
      auto dy_tlm = model.forward_tlm(x, eps*dx);

      std::cout << "\n=== TLM Finite Difference vs Autodiff Check ===\n";
      std::cout << "Finite Difference (FD):  " << dy_fd << "\n";
      std::cout << "Tangent Linear Model (TLM): " << dy_tlm << "\n";
      auto diff_tlm = torch::max(torch::abs(dy_fd - dy_tlm)).item<double>();
      std::cout << "Max Difference (FD vs TLM): " << diff_tlm << "\n";
      printTestResult("TLM finite difference test", diff_tlm < 1e-6);

      // Check linearity of TLM
      double scale = 10.5;
      auto y_tlm1 = model.forward_tlm(x, scale*dx);
      auto y_tlm2 = scale * model.forward_tlm(x, dx);
      auto diff_tlm2 = torch::max(torch::abs(y_tlm1 - y_tlm2)).item<double>();
      std::cout << "\n=== TLM Linearity Check ===\n";
      std::cout << "TLM1: " << y_tlm1 << "\n";
      std::cout << "TLM2: " << y_tlm2 << "\n";
      std::cout << "Max Difference (TLM1 vs TLM2): " << diff_tlm2 << "\n";
      printTestResult("TLM linearity test", diff_tlm2 < 1e-6);

      // ---------- AD Check: dot product test ----------
      auto tlm = model.forward_tlm(x, dx);
      auto adj = model.forward_ad(x, dy);

      double dot1 = dot_local(tlm, dy);
      double dot2 = dot_local(dx, adj);
      std::cout << "\n=== Adjoint Dot Product Test ===\n";
      std::cout << "<TLM dx, dy>: " << dot1 << "\n";
      std::cout << "<dx, AD dy>: " << dot2 << "\n";
      double relative_diff = std::abs(dot1 - dot2) / std::abs(dot1);
      printTestResult("Adjoint dot product test", relative_diff < 1e-6);

      return 0;
    }
   private:
    std::string appname() const {
      return "soca::test_OceanIceFFNN<";
    }
  };
}
