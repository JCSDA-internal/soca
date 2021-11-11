/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/AnalyticInit/AnalyticInit.h"
#include "soca/AnalyticInit/AnalyticInitFortran.h"

namespace soca {

  static oops::AnalyticInitMaker<ufo::ObsTraits, AnalyticInit>
    makerAnalyticInit_("soca_ana_init");

  void AnalyticInit::fillGeoVaLs(const ufo::Locations & locs,
                                 ufo::GeoVaLs & geovals) const {
    soca_analytic_geovals_f90(geovals.toFortran(), locs);
  }
}  // namespace soca
