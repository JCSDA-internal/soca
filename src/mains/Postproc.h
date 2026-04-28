/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "soca/Geometry/Geometry.h"
#include "soca/PostProcess/PostProcessIce.h"
#include "soca/State/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

namespace soca {

class Postproc : public oops::Application {
 public:
  // -----------------------------------------------------------------------------

  explicit Postproc(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
  static const std::string classname() {return "soca::Postproc";}

  // -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    const eckit::LocalConfiguration geomConfig(fullConfig, "geometry");
    const soca::Geometry geom(geomConfig, this->getComm());
    const eckit::LocalConfiguration restartConfig(fullConfig, "restart");
    const soca::State restart(geom, restartConfig);
    const eckit::LocalConfiguration analysisConfig(fullConfig, "analysis");
    const soca::State analysis(geom, analysisConfig);
    soca::State pproc(geom, restart);
    PostProcessIce ppIce(geom, fullConfig.getSubConfiguration("postprocess ice"));
    ppIce.postProcess(pproc, restart, analysis);
    if (fullConfig.has("restart output")) {
      const eckit::LocalConfiguration restartOutputConfig(fullConfig, "restart output");
      pproc.write(restartOutputConfig);
    }
    return 0;
  }

 private:
  // -----------------------------------------------------------------------------
  std::string appname() const {
    return "soca::Postproc";
  }
};

}  // namespace soca
