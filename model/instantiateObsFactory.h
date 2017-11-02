// Used in executable/MakeObs.cc

#ifndef MOM5CICE5_MODEL_INSTANTIATEOBSFACTORY_H_
#define MOM5CICE5_MODEL_INSTANTIATEOBSFACTORY_H_

#include "model/ObsOperator/ObsFraction.h"
//#include "model/ObsWindQG.h"
//#include "model/ObsWSpeedQG.h"
#include "model/ObsOperator/Observation.h"

namespace mom5cice5 {

  void instantiateObsFactory() {
    static ObsMaker<ObsFraction> makerFraction_("Fraction");
    //static ObsMaker<ObsWSpeedQG> makerWSpeed_("WSpeed");
  }

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_INSTANTIATEOBSFACTORY_H_
