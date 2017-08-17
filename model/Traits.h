
#ifndef MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_
#define MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_

#include <string>

#include "model/Geometry.h"
#include "model/Increment.h"
#include "model/State.h"
#include "model/Variables.h"
#include "model/ErrorCovariance.h"

namespace mom5cice5 {

struct Traits {
  static std::string name() {return "MOM5CICE5";}

  typedef mom5cice5::Geometry            Geometry;
  typedef mom5cice5::Variables           Variables;
  typedef mom5cice5::State               State;
  typedef mom5cice5::Increment           Increment;
  static std::string nameCovar() {return "MC5Error";}
  typedef mom5cice5::ErrorCovariance     Covariance;
  
};

}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_MOM5CICE5TRAITS_H_
