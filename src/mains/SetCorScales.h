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
      //  setup geometry
      const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
      const Geometry resol(resolConfig, this->getComm());

      //  rh from date and variables
      const util::DateTime thedate(fullConfig.getString("date"));
      const oops::Variables vars(fullConfig, "corr variables");
      Increment rh(resol, vars, thedate);

      //  compute horizontal decorrelation length scales
      const eckit::LocalConfiguration scalesConfig(fullConfig, "scales");
      rh.horiz_scales(scalesConfig);
      const eckit::LocalConfiguration rhoutputConfig(fullConfig, "rh output");
      rh.write(rhoutputConfig);
      oops::Log::test() << "Output horizontal scales: " << rh << std::endl;

      //  compute vertical decorrelation length scales
      const eckit::LocalConfiguration rvoutputConfig(fullConfig, "rv output");
      Increment rv(rh);
      const double vert = scalesConfig.getDouble("vert layers");
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
