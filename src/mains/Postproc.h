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
#include "soca/Increment/Increment.h"
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
  /// @brief Standalone CICE-restart postprocessing driver.
  ///
  /// Reads the background and the analysis increment, forms the analysis as
  /// `bg + incr`, then runs PostProcessIce on the per-category CICE restart.
  /// PostProcessIce owns the CICE restart read/write (paths from the
  /// `postprocess ice: cice restart:` block).
  int execute(const eckit::Configuration & fullConfig) const {
    const soca::Geometry geom(fullConfig.getSubConfiguration("geometry"),
                              this->getComm());

    std::unique_ptr<soca::State> analysis;
    if (fullConfig.has("analysis")) {
      analysis = std::make_unique<soca::State>(geom, fullConfig.getSubConfiguration("analysis"));
    } else {
      // Background -- aggregate ice + ocean (history-style ice file).
      const soca::State bg(geom, fullConfig.getSubConfiguration("background"));

      // 3DVar increment, read on the bg's variables and time.
      soca::Increment incr(geom, bg.variables(), bg.validTime());
      incr.read(fullConfig.getSubConfiguration("increment"));

      // analysis = bg + incr
      analysis = std::make_unique<soca::State>(bg);
      *analysis += incr;
    }
    // Per-category CICE restart, read by PostProcessIce itself. postProcess
    // internally copies ocean T/S off `analysis` onto pproc, so no caller-
    // side ocean merge is needed -- the SST-on-ice2noice adjustment writes
    // to pproc.T which starts at analysis.T.
    PostProcessIce ppIce(geom, fullConfig.getSubConfiguration("postprocess ice"));
    soca::State restart = ppIce.readRestart(geom, analysis->validTime());
    soca::State pproc(geom, restart);
    ppIce.postProcess(pproc, restart, *analysis);
    ppIce.writeRestart(pproc, analysis->validTime());
    return 0;
  }

 private:
  // -----------------------------------------------------------------------------
  std::string appname() const {
    return "soca::Postproc";
  }
};

}  // namespace soca
