/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_ANALYTICINIT_ANALYTICINIT_H_
#define SOCA_ANALYTICINIT_ANALYTICINIT_H_

#include "oops/interface/AnalyticInitBase.h"
#include "ufo/ObsTraits.h"

namespace soca {

  class AnalyticInitParameters : public oops::AnalyticInitParametersBase {
    OOPS_CONCRETE_PARAMETERS(AnalyticInitParameters, AnalyticInitParametersBase)
   public:
    // No parameters, for now
  };

// -----------------------------------------------------------------------------


  class AnalyticInit :
    public oops::interface::AnalyticInitBase<ufo::ObsTraits> {
   public:
    typedef AnalyticInitParameters Parameters_;

    explicit AnalyticInit(const Parameters_ &) {}
    void fillGeoVaLs(const ufo::Locations &, ufo::GeoVaLs &) const override;
  };

}  // namespace soca
#endif  // SOCA_ANALYTICINIT_ANALYTICINIT_H_
