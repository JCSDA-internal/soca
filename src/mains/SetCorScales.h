/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MAINS_SETCORSCALES_H_
#define MAINS_SETCORSCALES_H_

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/Point.h"

#include "soca/Traits.h"

#include "soca/Geometry/Geometry.h"
#include "soca/State/State.h"
#include "soca/Increment/Increment.h"

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/PostProcessor.h"
#include "oops/generic/gc99.h"
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

      // locate islands
      if (fullConfig.has("islands")) {
        oops::Log::info() << "====================== apply explicit island filtering" << std::endl;
        // Get the field set from rh
        atlas::FieldSet rhFs;
        rh.toFieldSet(rhFs);
        // Get the 2D grid
        auto lonlat = atlas::array::make_view<double, 2>(resol.functionSpace().lonlat());
        const eckit::LocalConfiguration islandsConfig(fullConfig, "islands");
        // Get the location of Islands from the config
        std::vector<atlas::PointLonLat> p0;
        const std::vector<double> lons = islandsConfig.getDoubleVector("lon");
        const std::vector<double> lats = islandsConfig.getDoubleVector("lat");
        ASSERT(lons.size() == lats.size());
        for (size_t i = 0; i < lons.size(); ++i) {
            p0.push_back(atlas::PointLonLat(lons[i], lats[i]));
        }
        // Recompute the correlation scales
        float scale = islandsConfig.getFloat("scale");
        for (auto point : p0) {
          oops::Log::info() << "---------- Island location: " << point << std::endl;
          for (auto & field : rhFs) {
            oops::Log::info() << "---------- Field name: " << field.name() << std::endl;
            auto view = atlas::array::make_view<double, 2>(field);
              for (int jnode = 0; jnode < field.shape(0); ++jnode) {
                atlas::PointLonLat p1(lonlat(jnode, 0), lonlat(jnode, 1));
                double d = atlas::util::Earth::distance(point, p1)/1000.0;
                for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
                  view(jnode, jlevel) *= (1.0 - oops::gc99(d/scale));
                }
              }
            }
          }
        }

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
