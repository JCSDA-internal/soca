/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
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
  /// `bg + incr`, then hands it to PostProcessIce::postprocess. The
  /// PostProcessIce call owns the CICE restart read/write (paths from the
  /// `postprocess ice: cice restart:` block) and returns an aggregate-ice
  /// State that this driver discards (the standalone binary only cares about
  /// the on-disk restart).
  int execute(const eckit::Configuration & fullConfig) const {
    const soca::Geometry geom(fullConfig.getSubConfiguration("geometry"),
                              this->getComm());

    std::unique_ptr<soca::State> analysis;
    if (fullConfig.has("analysis")) {
      analysis = std::make_unique<soca::State>(geom, fullConfig.getSubConfiguration("analysis"));
    } else {
      // Background -- aggregate ice + ocean (history-style ice file).
      const soca::State bg(geom, fullConfig.getSubConfiguration("background"));
      oops::Log::test() << " Background: " << bg << std::endl;

      // 3DVar increment, read on the bg's variables and time.
      soca::Increment incr(geom, bg.variables(), bg.validTime());
      incr.read(fullConfig.getSubConfiguration("increment"));

      // analysis = bg + incr
      analysis = std::make_unique<soca::State>(bg);
      *analysis += incr;
    }
    oops::Log::test() << " Analysis: " << *analysis << std::endl;
    PostProcessIce ppIce(geom, fullConfig.getSubConfiguration("postprocess ice"));
    const soca::State result = ppIce.postprocess(*analysis);
    oops::Log::test() << " Postprocessed state: " << result << std::endl;
    return 0;
  }

 private:
  // -----------------------------------------------------------------------------
  std::string appname() const {
    return "soca::Postproc";
  }
};

}  // namespace soca
