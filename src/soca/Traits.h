/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_TRAITS_H_
#define SOCA_TRAITS_H_

#include <string>

#include "soca/Geometry/Geometry.h"
#include "soca/Increment/Increment.h"
#include "soca/Covariance/ErrorCovariance.h"
#include "soca/GetValuesTraj/GetValuesTraj.h"
#include "soca/Localization/Localization.h"
#include "soca/ModelBias.h"
#include "soca/ModelBiasIncrement.h"
#include "soca/ModelBiasCovariance.h"
#include "soca/State/State.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"
#include "ufo/Locations.h"

namespace soca {

struct Traits {
  static std::string name() {return "SOCA";}
  static std::string nameCovar() {return "SocaError";}
  static std::string nameCovar4D() {return "SocaError";}

  typedef soca::Geometry            Geometry;
  typedef soca::State               State;
  typedef soca::Increment           Increment;
  typedef soca::ErrorCovariance     Covariance;
  typedef soca::GetValuesTraj       InterpolatorTraj;

  typedef soca::ModelBias           ModelAuxControl;
  typedef soca::ModelBiasIncrement  ModelAuxIncrement;
  typedef soca::ModelBiasCovariance ModelAuxCovariance;
  typedef soca::Localization        Localization;

  typedef ufo::ObsBias              ObsAuxControl;
  typedef ufo::ObsBiasCovariance    ObsAuxCovariance;
  typedef ufo::ObsBiasIncrement     ObsAuxIncrement;
  typedef ufo::ObsDiagnostics       ObsDiagnostics;
  //  typedef ufo::ObsCheck             ObsCheck;
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
