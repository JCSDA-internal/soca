/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MAINS_SETCORSCALES_H_
#define MAINS_SETCORSCALES_H_

#include <string>

#include "soca/Traits.h"

#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"
#include "soca/Increment/Increment.h"

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace soca {

  class SetCorScales : public oops::Application {
   public:
    explicit SetCorScales(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}
    static const std::string classname() {return "soca::SetCorScales";}

    int execute(const eckit::Configuration & fullConfig) const {
      //  Setup geometry
      const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
      const Geometry resol(resolConfig, this->getComm());

      //  Get date
      const util::DateTime thedate(fullConfig.getString("date"));

      //  Variables
      const oops::Variables vars(fullConfig, "corr variables");

      //  Get correlation scaling parameters
      const eckit::LocalConfiguration scalesConfig(fullConfig, "scales");
      const double r_mult = scalesConfig.getDouble("rossby mult");
      const double r_min_grid = scalesConfig.getDouble("min grid mult");
      const double vert = scalesConfig.getDouble("vert layers");

      //  Compute horizontal decorrelation length scales
      const eckit::LocalConfiguration rhoutputConfig(fullConfig, "rh output");
      Increment rh(resol, vars, thedate);
      rh.horiz_scales(r_mult, r_min_grid);
      rh.write(rhoutputConfig);
      oops::Log::test() << "Output horizontal scales: " << rh << std::endl;

      //  Compute vertical decorrelation length scales
      const eckit::LocalConfiguration rvoutputConfig(fullConfig, "rv output");
      Increment rv(rh);
      rv.vert_scales(vert);
      rv.write(rvoutputConfig);
      oops::Log::test() << "Output vertical scales: " << rv << std::endl;

      return 0;
    }
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "soca::SetCorScales<";
    }
    // -----------------------------------------------------------------------------
  };

}  // namespace soca
#endif  // MAINS_SETCORSCALES_H_
