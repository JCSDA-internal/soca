/*
 * (C) Copyright 2020 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Traits.h"
#include "ModelUFS.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace soca {
  class ModelBias;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

class ModelUFS: public oops::interface::ModelBase<Traits>,
                private util::ObjectCounter<ModelUFS> {
 public:
  static const std::string classname() {return "soca::ModelUFS";}

  ModelUFS(const Geometry &, const eckit::Configuration &);
  ~ModelUFS();

  void initialize(State &) const;
  void step(State &, const ModelBias &) const;

/// Finish model integration
  void finalize(State &) const;
  int saveTrajectory(State &, const ModelBias &) const;

  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry geom_;
  const oops::Variables vars_;
  char jedidir_[10000];
  char ufsdir_[10000];
};
// -----------------------------------------------------------------------------

}  // namespace soca
