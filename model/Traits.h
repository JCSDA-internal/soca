
#ifndef SOCA_MODEL_SOCATRAITS_H_
#define SOCA_MODEL_SOCATRAITS_H_

#include <string>

#include "model/ModelAtLocations/Gom.h"
#include "model/Locations/Loc.h"
#include "model/ModelBias.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelBiasCovariance.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsBiasCovariance.h"
#include "model/ObsSpace/ObsSpace.h"
#include "model/ObsVector/ObsVec.h"
//#include "model/Covariance/ErrorCovariance.h"
#include "model/Geometry/Geometry.h"
#include "model/Increment/Increment.h"
#include "model/Model/Model.h"
#include "model/State/State.h"
#include "model/Variables/Variables.h"

namespace soca {

struct Traits {
  static std::string name() {return "SOCA";}
  
  typedef soca::Geometry            Geometry;
  typedef soca::Variables           Variables;

  typedef soca::State               State;
  typedef soca::Model               Model;
  typedef soca::Increment           Increment;
  //static std::string nameCovar() {return "MC5Error";}
  //typedef soca::ErrorCovariance     Covariance;

  typedef soca::ModelBias           ModelAuxControl;
  typedef soca::ModelBiasIncrement  ModelAuxIncrement;
  typedef soca::ModelBiasCovariance ModelAuxCovariance;
  
  typedef soca::ObsSpace            ObsSpace;
  //typedef soca::Observation         ObsOperator;
  //  typedef soca::LinearObsOp         LinearObsOperator;
  typedef soca::ObsVec              ObsVector;

  typedef soca::ObsBias             ObsAuxControl;
  typedef soca::ObsBiasIncrement    ObsAuxIncrement;
  typedef soca::ObsBiasCovariance   ObsAuxCovariance;
  
  typedef soca::Gom                 GeoVaLs; //ModelAtLocations;
  typedef soca::Loc                 Locations;
};

}  // namespace soca

#endif  // SOCA_MODEL_SOCATRAITS_H_
