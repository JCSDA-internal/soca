/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"

#include "soca/MLBalance/KEmul/IceEmul.h"

namespace soca {
  class IceEmulApp : public oops::Application {
   public:
    explicit IceEmulApp(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}

    // -----------------------------------------------------------------------------
    static const std::string classname() {return "soca::IceEmulApp";}

    // -----------------------------------------------------------------------------
    int execute(const eckit::Configuration & config) const {
      oops::Log::info() << "Initialize the FFNN" << std::endl;
      soca::IceEmul iceEmul(config, getComm());

      // Generate patterns-targets pairs and train
      if (config.has("training")) {
        oops::Log::info() << "Prepare patterns/targets pairs" << std::endl;
        std::string fileName;
        config.get("training.cice history", fileName);
        auto result = iceEmul.prepData(fileName);
        torch::Tensor inputs = std::get<0>(result);
        torch::Tensor targets = std::get<1>(result);

        oops::Log::info() << "Initialize the normalization" << std::endl;
        torch::Tensor mean = std::get<4>(result);
        torch::Tensor std = std::get<5>(result);
        iceEmul.getModel()->initNorm(mean, std);

        oops::Log::info() << "Train the FFNN" << std::endl;
        iceEmul.train(inputs, targets);
      }

      // Predictions
      if (config.has("prediction")) {
        // Assert that this is a serial job
        if (getComm().size() > 1) {
          oops::Log::error() << "Prediction is only supported in serial mode" << std::endl;
          return 1;
        }
        oops::Log::info() << "Predict" << std::endl;
        std::string fileName;
        config.get("prediction.cice history", fileName);
        std::string fileNameResults;
        config.get("prediction.output filename", fileNameResults);
        iceEmul.predict(fileName, fileNameResults);
      }
      return 0;
    }
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "soca::iceEmul";
    }
  };
}  // namespace soca

int main(int argc, char* argv[]) {
  oops::Run run(argc, argv);
  soca::IceEmulApp iceemulapp;
  return run.execute(iceemulapp);
}
