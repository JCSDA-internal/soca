
#ifndef MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_
#define MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_

#include <string>

#include "model/Gom.h"
#include "model/LinearObsOp.h"
#include "model/Loc.h"
#include "model/ObsSpace.h"
#include "model/ObsVec.h"
#include "model/Geometry.h"
#include "model/Increment.h"
#include "model/Model.h"
#include "model/ModelBias.h"
#include "model/Observation.h"
#include "model/State.h"
#include "model/Variables.h"
#include "model/ErrorCovariance.h"

namespace mom5cice5 {

struct Traits {
  static std::string name() {return "MOM5CICE5";}

  typedef mom5cice5::Geometry            Geometry;
  typedef mom5cice5::Variables           Variables;

  typedef mom5cice5::State               State;
  typedef mom5cice5::Model               Model;
  typedef mom5cice5::Increment           Increment;
  static std::string nameCovar() {return "MC5Error";}
  typedef mom5cice5::ErrorCovariance     Covariance;

  typedef mom5cice5::ModelBias           ModelAuxControl;
  
  typedef mom5cice5::ObsSpace            ObsSpace;
  typedef mom5cice5::Observation         ObsOperator;
  typedef mom5cice5::LinearObsOp         LinearObsOperator;
  typedef mom5cice5::ObsVec              ObsVector;

  typedef mom5cice5::Gom                 ModelAtLocations;
  typedef mom5cice5::Loc                 Locations;
};

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_
