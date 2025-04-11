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
  if (ii <= 0 || ii > static_cast<int>(traj_.size())) {
    throw eckit::OutOfRange("Index out of bounds in ModelTrajectory::get", Here());
  }
  return traj_[ii-1];
}
// -----------------------------------------------------------------------------

}  // namespace soca

