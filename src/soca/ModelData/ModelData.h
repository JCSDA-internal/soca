/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>
#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace soca {
  class Geometry;
}

// -----------------------------------------------------------------------------

namespace soca {

class ModelData : public util::Printable {
 public:
  static const std::string classname() {return "soca::ModelData";}
  oops::Variables defaultVariables() const {
    return oops::Variables(std::vector<std::string>({"skin_temperature_at_surface_where_sea",
            "sea_water_potential_temperature", "sea_water_salinity", "sea_water_cell_thickness",
            "distance_from_coast", "sea_area_fraction", "sea_ice_snow_thickness",
            "sea_surface_height_above_geoid", "sea_ice_area_fraction",
            "sea_ice_thickness", "Carbon_nitrogen_detritus_concentration",
            "Particulate_inorganic_carbon", "colored_dissolved_organic_carbon",
            "diatom_concentration", "chlorophyte_concentration", "cyano-bacteria_concentration",
            "coccolithophore_concentration", "dinoflagellate_concentration",
            "phaeocystis_concentration"}));
  }

  explicit ModelData(const Geometry &) {}
  ~ModelData() = default;

  const eckit::LocalConfiguration modelData() const {return eckit::LocalConfiguration();}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace soca
