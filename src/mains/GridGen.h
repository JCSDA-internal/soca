/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MAINS_GRIDGEN_H_
#define MAINS_GRIDGEN_H_

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/PostProcessor.h"
#include "soca/Geometry/Geometry.h"
#include "soca/Model/Model.h"
#include "oops/runs/Application.h"


namespace soca {

  class GridGen : public oops::Application {
   public:
    static const std::string classname() {return "soca::GridGen";}

    int execute(const eckit::Configuration & fullConfig) const {
      //  Setup resolution
      const eckit::LocalConfiguration geomconfig(fullConfig, "geometry");
      const Geometry geom(geomconfig);

      //  Generate model grid
      const eckit::LocalConfiguration gridgenconfig(fullConfig, "gridgen");
      geom.gridgen(gridgenconfig);

      return 0;
    }
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "soca::GridGen<";
    }
    // -----------------------------------------------------------------------------
  };

}  // namespace soca
#endif  // MAINS_GRIDGEN_H_
