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

        // Read the backgrounds
        oops::Log::info() << "====================== read bkg" << std::endl;
        // Loop through the list of backgrounds

        const eckit::LocalConfiguration bkgConfig(config, "backgrounds");
        const std::vector<eckit::LocalConfiguration> bkgList = bkgConfig.getSubConfigurations();
        oops::Log::info() << "Number of backgrounds: " << bkgList.size() << std::endl;
        for (const eckit::LocalConfiguration & bkgConf : bkgList) {
          soca::State xb(geom, config.socaVars, config.cycleDate);
          xb.read(bkgConf);
          atlas::FieldSet xbFs;
          xb.toFieldSet(xbFs);
          oops::Log::info() << "Background:" << std::endl;
          oops::Log::info() << xb << std::endl;
        }

        // Prepare the patterns/targets pairs
        oops::Log::info() << "====================== prepare patterns/targets" << std::endl;
        // Open the torch output file
        const std::string outputFile = config.getString("output file");
        torch::Tensor patterns, targets;

        // Assuming patterns and targets are filled with the appropriate data
        // This is just a placeholder for the actual data preparation logic
        patterns = torch::rand({10, 3, 224, 224}); // Example tensor
        targets = torch::rand({10, 1}); // Example tensor

        // Save the tensors to a file
        torch::save({patterns, targets}, outputFile);
        oops::Log::info() << "Patterns and targets saved to " << outputFile << std::endl;

       return 0;
      }
  };
}