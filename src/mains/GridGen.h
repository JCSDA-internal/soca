/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MAINS_GRIDGEN_H_
#define MAINS_GRIDGEN_H_

#include <string>

#include "soca/Traits.h"

#include "soca/Geometry/Geometry.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/base/PostProcessor.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

namespace soca {

  class GridGen : public oops::Application {
   public:
    explicit GridGen(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}
    static const std::string classname() {return "soca::GridGen";}

    int execute(const eckit::Configuration & fullConfig) const {
      //  Setup resolution
      const eckit::LocalConfiguration geomconfig(fullConfig, "geometry");
      const Geometry geom(geomconfig, this->getComm(), true);
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
