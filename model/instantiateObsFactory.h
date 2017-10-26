
#ifndef MOM5CICE5_MODEL_INSTANTIATEOBSFACTORY_H_
#define MOM5CICE5_MODEL_INSTANTIATEOBSFACTORY_H_

#include "model/ObsFraction.h"
//#include "model/ObsWindQG.h"
//#include "model/ObsWSpeedQG.h"
#include "model/Observation.h"

namespace mom5cice5 {

  void instantiateObsFactory() {
    static ObsMaker<ObsFraction> makerFraction_("Fraction");
    //static ObsMaker<ObsWSpeedQG> makerWSpeed_("WSpeed");
  }

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_INSTANTIATEOBSFACTORY_H_
