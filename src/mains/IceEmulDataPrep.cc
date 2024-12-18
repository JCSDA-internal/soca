/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <torch/torch.h>
#include <torch/csrc/distributed/c10d/ProcessGroupMPI.hpp>
#include <torch/csrc/distributed/c10d/Utils.hpp>

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"

#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"

namespace soca {
  class IceEmulDataPrepApp : public oops::Application {
   public:
    explicit IceEmulDataPrepApp(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}

    // -----------------------------------------------------------------------------
    static const std::string classname() {return "soca::IceEmulDataPrepApp";}

    // -----------------------------------------------------------------------------
    int execute(const eckit::Configuration & config, bool /*validate*/) const {
      oops::Log::info() << "Read backgrounds and prepare the patterns/targets pairs" << std::endl;

      // Setup the soca geometry
      oops::Log::info() << "====================== geometry" << std::endl;
      const eckit::LocalConfiguration geomConfig(config, "geometry");
      const soca::Geometry geom(geomConfig, this->getComm());

      // Get the list of input and output variables from the configuration
      const std::vector<std::string> inputVars = config.getStringVector("input variables");
      const std::vector<std::string> outputVars = config.getStringVector("output variables");
      const int outputSize = outputVars.size();
      const int inputSize = inputVars.size();

      // Read the backgrounds
      oops::Log::info() << "====================== read bkg" << std::endl;
      const eckit::LocalConfiguration bkgConfig(config, "backgrounds");
      const std::vector<eckit::LocalConfiguration> bkgList = bkgConfig.getSubConfigurations();

      // Temporary containers for patterns/targets pairs. No real equivalent to push_back for torch tensors
      std::vector<std::vector<float>> patternsVec(inputSize);
      std::vector<std::vector<float>> targetsVec(outputSize);

      // Loop through the backgrounds
      oops::Log::info() << "Number of backgrounds: " << bkgList.size() << std::endl;
      for (const eckit::LocalConfiguration & bkgConf : bkgList) {
        soca::State xb(geom, bkgConf);
        atlas::FieldSet xbFs;
        xb.toFieldSet(xbFs);
        oops::Log::info() << "Background:" << std::endl;
        oops::Log::info() << xb << std::endl;

        // Get the sea ice area fraction, used to filter the patterns
        auto viewAice = atlas::array::make_view<double, 2>(xbFs["sea_ice_area_fraction"]);

        // Loop through the input variables
        int varIndex = 0;
        for ( auto inputVar : inputVars ) {
          auto view = atlas::array::make_view<double, 2>(xbFs[inputVar]);

          // Store the patterns only where the sea ice area fraction is between 0 and 1
          for (int i = 0; i < view.shape(0); ++i) {
            for (int j = 0; j < view.shape(1); ++j) {
              if (viewAice(i, j) > 0.0 && viewAice(i, j) <= 1.0) {
                patternsVec[varIndex].push_back(view(i, j));
              }
            }
          }
          varIndex++;
        }  // end loop through input variables

        // Loop through the output variables
        varIndex = 0;
        for ( auto outputVar : outputVars ) {
          auto view = atlas::array::make_view<double, 2>(xbFs[outputVar]);

          // Store the targets only where the sea ice area fraction is between 0 and 1
          for (int i = 0; i < view.shape(0); ++i) {
            for (int j = 0; j < view.shape(1); ++j) {
              if (viewAice(i, j) > 0.0 && viewAice(i, j) <= 1.0) {
                targetsVec[varIndex].push_back(view(i, j));
              }
            }
          }
          varIndex++;
        }  // end loop through output variables
      }  // end loop through backgrounds

      // Prepare the patterns/targets pairs
      oops::Log::info() << "====================== prepare patterns/targets" << std::endl;

      const int numPatterns = patternsVec.empty() ? 0 : patternsVec[1].size();
      std::cout << "patternsVec size: " << numPatterns << std::endl;
      torch::Tensor patterns = torch::empty({numPatterns, inputSize});
      torch::Tensor targets = torch::empty({numPatterns, outputSize});

      // Fill the patterns and targets tensors
      for (int i = 0; i < numPatterns; ++i) {
        for (int j = 0; j < inputSize; ++j) {
          patterns[i][j] = patternsVec[j][i];
        }
        for (int j = 0; j < outputSize; ++j) {
          targets[i][j] = patternsVec[j][i];
        }
      }

      // Initialize torch distributed
      //torch::distributed::c10::init_process_group(torch::distributed::Backend::MPI);
      torch::distributed::c10d::ProcessGroupMPI::Options options;
      //auto process_group = torch::distributed::c10d::ProcessGroupMPI::createProcessGroupMPI(options);

      /*
      // Gather all patterns and targets tensors from all ranks
      std::vector<torch::Tensor> gatheredPatterns;
      std::vector<torch::Tensor> gatheredTargets;
      torch::distributed::c10::all_gather(gatheredPatterns, patterns);
      torch::distributed::c10::all_gather(gatheredTargets, targets);

      // Concatenate the vector of tensors into a single tensor
      torch::Tensor allPatterns = torch::cat(gatheredPatterns, 0);
      torch::Tensor allTargets = torch::cat(gatheredTargets, 0);

      // Save the concatenated tensors to file only on rank 0
      if (this->getComm().rank() == 0) {
        const std::string outputFile = config.getString("output file");
        torch::save({allPatterns, allTargets}, outputFile);
        oops::Log::info() << "Patterns and targets saved to " << outputFile << std::endl;
      }
      */
      return 0;
    }
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "soca::IceEmulDataPrepApp";
    }
  };
}  // namespace soca

int main(int argc, char* argv[]) {
  oops::Run run(argc, argv);
  soca::IceEmulDataPrepApp app;
  return run.execute(app);
}
