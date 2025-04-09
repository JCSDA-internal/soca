/*
* (C) Copyright 2025 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include "soca/LinearModel/OceanIceEmulator/ModelTrajectory.h"

#include "eckit/exception/Exceptions.h"
#include "soca/State/State.h"

// -----------------------------------------------------------------------------
namespace soca {
// -----------------------------------------------------------------------------
ModelTrajectory::ModelTrajectory(const bool ltraj) : ltraj_(ltraj), traj_() {}
// -----------------------------------------------------------------------------
ModelTrajectory::~ModelTrajectory() {}
// -----------------------------------------------------------------------------
void ModelTrajectory::set(const State & xx) {
  if (ltraj_) traj_.push_back(new State(xx));
}
// -----------------------------------------------------------------------------
const State & ModelTrajectory::get(const int ii) const {
  // TODO(G): Pass the max number of elements in the trajectory as a parameter
  //if (ii < 1 || ii > 2) {
  //    throw eckit::BadValue("ModelTrajectory::get index out of range", Here());
  //}
  //ASSERT(ltraj_);
  //ASSERT(traj_.size() == 2);
  //ASSERT(1 <= ii && ii <= 2);
  return traj_[ii-1];
}
// -----------------------------------------------------------------------------

}  // namespace soca

