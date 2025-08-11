/*
* (C) Copyright 2025 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#ifndef SOCA_LINEARMODEL_OCEANICEEMULATOR_LINEARMODELOCEANICEEMULATOR_H_
#define SOCA_LINEARMODEL_OCEANICEEMULATOR_LINEARMODELOCEANICEEMULATOR_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "soca/Geometry/Geometry.h"
#include "soca/Model/OceanIceEmulator/ModelOceanIceEmulator.h"
#include "soca/LinearModel/OceanIceEmulator/ModelTrajectory.h"
#include "soca/ModelBias/ModelBias.h"
#include "soca/State/State.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace soca {
  class Geometry;
  class ModelBias;
  class ModelOceanIceEmulator;
  class ModelTrajectory;
  class State;
  class Increment;
}

// -----------------------------------------------------------------------------

namespace soca {

// -----------------------------------------------------------------------------
/// OceanIceEmulator linear model definition.
/*!
 *  OceanIceEmulator linear model definition and configuration parameters.
 */

class LinearModelOceanIceEmulator: public util::Printable,
             private util::ObjectCounter<LinearModelOceanIceEmulator>
{
 public:
  static const std::string classname() {return "soca::LinearModelOceanIceEmulator";}
  static std::vector<std::string> names() {return {"LinearModelOceanIceEmulator"};}

  LinearModelOceanIceEmulator(const Geometry &, const eckit::Configuration &);
  ~LinearModelOceanIceEmulator();

  /// Prepare model integration
  void initializeTL(Increment &) const;
  void initializeAD(Increment &) const;

  /// Model integration
  void stepTL(Increment &, const ModelBiasIncrement &) const;
  void stepAD(Increment &, ModelBiasIncrement &) const;
  void setTrajectory(const State &, State &, const ModelBias &);

  /// Finish model integration
  void finalizeTL(Increment &) const;
  void finalizeAD(Increment &) const;

  /// Other utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return tstep_;}
  const Geometry & resolution() const {return geom_;}

 private:
  const ModelTrajectory * getTrajectory(const util::DateTime &) const;
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, ModelTrajectory * >::const_iterator trajICst;

  // Data
  util::Duration tstep_;
  const Geometry & geom_;
  const util::Duration steptraj_;
  std::map< util::DateTime, ModelTrajectory * > traj_;
  const soca::ModelOceanIceEmulator lrmodel_;
  const oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_LINEARMODEL_OCEANICEEMULATOR_LINEARMODELOCEANICEEMULATOR_H_
