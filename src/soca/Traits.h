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
// Geometry key type
typedef int F90geom;
// Geometry iterator key type
typedef int F90iter;
// Model key type
typedef int F90model;
// Tlm key type
typedef int F90tlm;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
typedef int F90bmat;
// ObOp trajectory
typedef int F90ootrj;
// State key
typedef int F90state;
// Increment key
typedef int F90inc;
// Variable change
typedef int F90varcha;
// GetValues key
typedef int F90getvalues;
typedef int F90lineargetvalues;


}  // namespace soca

#endif  // SOCA_TRAITS_H_
