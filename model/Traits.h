/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_SOCATRAITS_H_
#define SOCA_MODEL_SOCATRAITS_H_

#include <string>

#include "model/ModelBias.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelBiasCovariance.h"
#include "model/Geometry/Geometry.h"
#include "model/Increment/Increment.h"
#include "model/LocalizationMatrix/LocalizationMatrix.h"
#include "model/Model/Model.h"
#include "model/State/State.h"
#include "model/Covariance/ErrorCovariance.h"
#include "model/ObsBiasCovariance.h"
#include "model/Nothing/Nothing.h"

#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsCheck.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsOperator.h"
#include "ufo/LinearObsOperator.h"

namespace soca {

struct Traits {
  static std::string name() {return "SOCA";}
  static std::string nameCovar() {return "SocaError";}
  
  typedef soca::Geometry            Geometry;
  typedef soca::State               State;
  typedef soca::Model               Model;
  typedef soca::Increment           Increment;
  typedef soca::ErrorCovariance     Covariance;
  typedef soca::Nothing               InterpolatorTraj;
  
  typedef soca::ModelBias           ModelAuxControl;
  typedef soca::ModelBiasIncrement  ModelAuxIncrement;
  typedef soca::ModelBiasCovariance ModelAuxCovariance;
  typedef soca::LocalizationMatrix  LocalizationMatrix;  

  //typedef soca::ObsBiasCovariance     ObsAuxCovariance;
  
  typedef ufo::ObsBias              ObsAuxControl;
  typedef ufo::ObsBiasIncrement     ObsAuxIncrement;
  typedef ufo::ObsBiasCovariance    ObsAuxCovariance;
  typedef ufo::ObsCheck             ObsCheck;
  typedef ufo::ObsOperator          ObsOperator;
  typedef ufo::LinearObsOperator    LinearObsOperator;
  typedef ufo::GeoVaLs              GeoVaLs;
  typedef ioda::Locations           Locations;
  typedef ioda::ObsSpace             ObsSpace;
  typedef ioda::ObsVector            ObsVector;
};

}  // namespace soca

#endif  // SOCA_MODEL_SOCATRAITS_H_
