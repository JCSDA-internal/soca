/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_ANALYTICINIT_ANALYTICINIT_H_
#define SOCA_ANALYTICINIT_ANALYTICINIT_H_

#include <string>

#include "oops/interface/AnalyticInitBase.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/ObsTraits.h"

namespace soca {

  class AnalyticInitParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(AnalyticInitParameters, Parameters)
    oops::RequiredParameter<std::string> method{"method", this};
   public:
    // No parameters, for now
  };

// -----------------------------------------------------------------------------


  class AnalyticInit :
    public oops::interface::AnalyticInitBase<ufo::ObsTraits> {
   public:
    explicit AnalyticInit(const eckit::Configuration &) {}
    void fillGeoVaLs(const ufo::SampledLocations &, ufo::GeoVaLs &) const override;
  };

}  // namespace soca
#endif  // SOCA_ANALYTICINIT_ANALYTICINIT_H_
