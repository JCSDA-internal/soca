/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MAINS_CHECKPOINTMODEL_H_
#define MAINS_CHECKPOINTMODEL_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/PostProcessor.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Model/Model.h"
#include "soca/State/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace soca {

  class CheckpointModel : public oops::Application {
   public:
    static const std::string classname() {return "soca::CheckpointModel";}

    int execute(const eckit::Configuration & fullConfig) const {
      //  Setup resolution
      const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
      const Geometry resol(resolConfig);

      //  Setup Model
      const eckit::LocalConfiguration modelConfig(fullConfig, "model");
      const Model model(resol, modelConfig);

      //  Setup state to write in the restart
      const eckit::LocalConfiguration backgroundConfig(fullConfig,
                                                       "background");
      State xb(resol, model.variables(), backgroundConfig);
      oops::Log::test() << "input background: " << std::endl << xb << std::endl;

      //  Setup state to write in the restart
      const eckit::LocalConfiguration analysisConfig(fullConfig, "analysis");
      State xa(resol, model.variables(), analysisConfig);
      oops::Log::test() << "analysis: " << std::endl << xa << std::endl;

      //  Initialize model
      model.initialize(xb);

      //  Set background to analysis
      xb = xa;

      //  Finalize model (dump restart)
      model.finalize(xb);

      oops::Log::test() << "output background: " << std::endl << xb
                        << std::endl;

      return 0;
    }
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "soca::CheckpointModel<";
    }
    // -----------------------------------------------------------------------------
  };

}  // namespace soca
#endif  // MAINS_CHECKPOINTMODEL_H_
