
#ifndef SOCA_MODEL_SOCATRAITS_H_
#define SOCA_MODEL_SOCATRAITS_H_

#include <string>

//#include "model/ModelAtLocations/Gom.h"
//#include "model/Locations/Loc.h"
#include "model/ModelBias.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelBiasCovariance.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsBiasCovariance.h"
//#include "model/ObsSpace/ObsSpace.h"
//#include "model/ObsVector/ObsVec.h"
//#include "model/Covariance/ErrorCovariance.h"
#include "model/Geometry/Geometry.h"
#include "model/Increment/Increment.h"
#include "model/Model/Model.h"
#include "model/State/State.h"
#include "model/Variables/Variables.h"
#include "model/Run/Run.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsSpace.h"
#include "ufo/ObsVector.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsSpace;
  class ObsVector;
  class ObsBias;
  class ObsBiasIncrement;  
}


namespace soca {

struct Traits {
  static std::string name() {return "SOCA";}
  
  typedef soca::Geometry            Geometry;
  //typedef soca::Variables           Variables;

  typedef soca::Run                 Run;

  typedef soca::State               State;
  typedef soca::Model               Model;
  typedef soca::Increment           Increment;
  //static std::string nameCovar() {return "MC5Error";}
  //typedef soca::ErrorCovariance     Covariance;

  typedef soca::ModelBias           ModelAuxControl;
  typedef soca::ModelBiasIncrement  ModelAuxIncrement;
  typedef soca::ModelBiasCovariance ModelAuxCovariance;
  
  //typedef soca::Observation         ObsOperator;
  //  typedef soca::LinearObsOp         LinearObsOperator;

  typedef ufo::ObsBias             ObsAuxControl;
  typedef ufo::ObsBiasIncrement    ObsAuxIncrement;
  typedef soca::ObsBiasCovariance   ObsAuxCovariance;
  
  typedef ufo::GeoVaLs              GeoVaLs; //ModelAtLocations;
  typedef ufo::Locations            Locations;
  typedef ufo::ObsSpace             ObsSpace;
  typedef ufo::ObsVector            ObsVector;
  
};

}  // namespace soca

#endif  // SOCA_MODEL_SOCATRAITS_H_
