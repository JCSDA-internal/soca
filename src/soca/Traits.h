/*
 * (C) Copyright 2017-2021 UCAR
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
#include "soca/LinearVariableChange/LinearVariableChange.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/ModelBias/ModelBiasCovariance.h"
#include "soca/ModelBias/ModelBiasIncrement.h"
#include "soca/State/State.h"
#include "soca/VariableChange/VariableChange.h"

namespace soca {

/**
 * \brief The main traits structure for SOCA.
 *
 * This structure is responsible for supplying SOCA specific code to the JEDI
 * applications within \ref src/mains directory.
 */
struct Traits {
  static std::string name() {return "SOCA";}
  static std::string nameCovar() {return "SocaError";}

  typedef soca::Geometry             Geometry;
  typedef soca::GeometryIterator     GeometryIterator;
  typedef soca::State                State;
  typedef soca::Increment            Increment;
  typedef soca::ErrorCovariance      Covariance;
  typedef soca::GetValues            GetValues;
  typedef soca::LinearGetValues      LinearGetValues;

  typedef soca::ModelBias            ModelAuxControl;
  typedef soca::ModelBiasIncrement   ModelAuxIncrement;
  typedef soca::ModelBiasCovariance  ModelAuxCovariance;

  typedef soca::LinearVariableChange LinearVariableChange;
  typedef soca::VariableChange       VariableChange;

};

}  // namespace soca

#endif  // SOCA_TRAITS_H_
