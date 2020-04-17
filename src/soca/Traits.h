/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRAITS_H_
#define SOCA_TRAITS_H_

#include <string>

#include "soca/Covariance/ErrorCovariance.h"
#include "soca/Geometry/Geometry.h"
#include "soca/GeometryIterator/GeometryIterator.h"
#include "soca/GetValues/GetValues.h"
#include "soca/GetValues/LinearGetValues.h"
#include "soca/Increment/Increment.h"
#include "soca/Localization/Localization.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/ModelBias/ModelBiasCovariance.h"
#include "soca/ModelBias/ModelBiasIncrement.h"
#include "soca/State/State.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"

namespace soca {

struct Traits {
  static std::string name() {return "SOCA";}
  static std::string nameCovar() {return "SocaError";}
  static std::string nameCovar4D() {return "SocaError";}

  typedef soca::Geometry            Geometry;
  typedef soca::GeometryIterator    GeometryIterator;
  typedef soca::State               State;
  typedef soca::Increment           Increment;
  typedef soca::ErrorCovariance     Covariance;
  typedef soca::GetValues           GetValues;
  typedef soca::LinearGetValues     LinearGetValues;

  typedef soca::ModelBias           ModelAuxControl;
  typedef soca::ModelBiasIncrement  ModelAuxIncrement;
  typedef soca::ModelBiasCovariance ModelAuxCovariance;
  typedef soca::Localization        Localization;

  typedef ufo::ObsBias              ObsAuxControl;
  typedef ufo::ObsBiasCovariance    ObsAuxCovariance;
  typedef ufo::ObsBiasIncrement     ObsAuxIncrement;
  typedef ufo::ObsDiagnostics       ObsDiagnostics;
  typedef ufo::ObsOperator          ObsOperator;
  typedef ufo::LinearObsOperator    LinearObsOperator;
  typedef ufo::GeoVaLs              GeoVaLs;
  typedef ufo::Locations            Locations;
  typedef ioda::ObsSpace            ObsSpace;
  typedef ioda::ObsVector           ObsVector;
  template <typename DATA> using ObsDataVector = ioda::ObsDataVector<DATA>;
};

}  // namespace soca

#endif  // SOCA_TRAITS_H_
