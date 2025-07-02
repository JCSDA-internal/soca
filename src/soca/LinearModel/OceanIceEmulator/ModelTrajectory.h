/*
* (C) Copyright 2025 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once


#include <string>

#include <boost/ptr_container/ptr_vector.hpp>

#include "oops/util/ObjectCounter.h"

namespace soca {
  class State;

/// OceanIceEmulator model trajectory

// -----------------------------------------------------------------------------
class ModelTrajectory: private util::ObjectCounter<ModelTrajectory> {
 public:
  static const std::string classname() {return "soca::ModelTrajectory";}

/// Constructor, destructor
  explicit ModelTrajectory(const bool ltraj = true);
  ~ModelTrajectory();

/// Save trajectory
  void set(const State &);

/// Get trajectory
  const State & get(const int) const;

 private:
  const bool ltraj_;
  boost::ptr_vector<State> traj_;
};
// -----------------------------------------------------------------------------

}  // namespace soca

