/*
* (C) Copyright 2025 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#ifndef SOCA_MODEL_OCEANICEEMULATOR_MODELOCEANICEEMULATOR_H_
#define SOCA_MODEL_OCEANICEEMULATOR_MODELOCEANICEEMULATOR_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}
namespace soca {
  class Geometry;
  class ModelBias;
  class State;
}

// -----------------------------------------------------------------------------

namespace soca {

  /// SOCA model definition.
  /*!
   *  SOCA nonlinear interface model definition and configuration parameters.
   */

  class ModelOceanIceEmulator : public util::Printable,
              private util::ObjectCounter<ModelOceanIceEmulator>
  {
   public:
    static const std::string classname() {return "soca::ModelOceanIceEmulator";}
    static std::vector<std::string> names() {return {"ModelOceanIceEmulator"};}

    ModelOceanIceEmulator(const Geometry &, const eckit::Configuration &);
    ~ModelOceanIceEmulator();

    /// Prepare model integration
    void initialize(State &) const;

    /// Model integration
    void step(State &, const ModelBias &) const;

    /// Finish model integration
    void finalize(State &) const;

    /// Utilities
    const util::Duration & timeResolution() const {return tstep_;}
    const oops::Variables & variables() const {return vars_;}

   private:
    void print(std::ostream &) const override;
    int keyConfig_;
    util::Duration tstep_;
    const Geometry & geom_;
    const oops::Variables vars_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_MODEL_OCEANICEEMULATOR_MODELOCEANICEEMULATOR_H_
