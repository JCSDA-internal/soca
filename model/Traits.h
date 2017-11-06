
#ifndef MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_
#define MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_

#include <string>

#include "model/ModelAtLocations/Gom.h"
#include "model/ObsOperator/LinearObsOp.h"
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
#include "model/ObsOperator/Observation.h"
#include "model/State/State.h"
#include "model/Variables/Variables.h"

namespace mom5cice5 {

struct Traits {
  static std::string name() {return "MOM5CICE5";}

  typedef mom5cice5::Geometry            Geometry;
  typedef mom5cice5::Variables           Variables;

  typedef mom5cice5::State               State;
  typedef mom5cice5::Model               Model;
  typedef mom5cice5::Increment           Increment;
  //static std::string nameCovar() {return "MC5Error";}
  //typedef mom5cice5::ErrorCovariance     Covariance;

  typedef mom5cice5::ModelBias           ModelAuxControl;
  typedef mom5cice5::ModelBiasIncrement  ModelAuxIncrement;
  typedef mom5cice5::ModelBiasCovariance ModelAuxCovariance;
  
  typedef mom5cice5::ObsSpace            ObsSpace;
  typedef mom5cice5::Observation         ObsOperator;
  typedef mom5cice5::LinearObsOp         LinearObsOperator;
  typedef mom5cice5::ObsVec              ObsVector;

  typedef mom5cice5::ObsBias             ObsAuxControl;
  typedef mom5cice5::ObsBiasIncrement    ObsAuxIncrement;
  typedef mom5cice5::ObsBiasCovariance   ObsAuxCovariance;
  
  typedef mom5cice5::Gom                 ModelAtLocations;
  typedef mom5cice5::Loc                 Locations;
};

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_
