/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <string>
#include <vector>

#include "saber/blocks/SaberBlockParametersBase.h"

namespace soca {

// --------------------------------------------------------------------------------------

class MLBalanceParameters : public saber::SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(MLBalanceParameters, saber::SaberBlockParametersBase)
 public:
  oops::RequiredParameter<eckit::LocalConfiguration> mlbalances{"ML Balances", this};
  oops::Variables mandatoryActiveVars() const override {
    return oops::Variables({"sea_water_potential_temperature",
                             "sea_water_salinity",
                             "sea_surface_height_above_geoid",
                             "sea_ice_area_fraction",
                             "sea_ice_thickness",
                             "sea_ice_snow_thickness"});}
};

// --------------------------------------------------------------------------------------
}  // namespace soca
