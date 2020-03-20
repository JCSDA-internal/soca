/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_MODEL_H_
#define SOCA_MODEL_MODEL_H_

#include <memory>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "soca/Fortran.h"

#include "oops/base/ModelBase.h"
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
  struct Traits;
}

// -----------------------------------------------------------------------------

namespace soca {

  /// SOCA model definition.
  /*!
   *  SOCA nonlinear model definition and configuration parameters.
   */

  class Model:public oops::ModelBase<Traits>, private util::ObjectCounter<Model>
  {
   public:
    static const std::string classname() {return "soca::Model";}

    Model(const Geometry &, const eckit::Configuration &);
    ~Model();

    /// Prepare model integration
    void initialize(State &) const;

    /// Model integration
    void step(State &, const ModelBias &) const;
    int saveTrajectory(State &, const ModelBias &) const;

    /// Finish model integration
    void finalize(State &) const;

    /// Utilities
    const util::Duration & timeResolution() const {return tstep_;}
    const oops::Variables & variables() const {return vars_;}

   private:
    void print(std::ostream &) const;
    int keyConfig_;
    util::Duration tstep_;
    bool setup_mom6_;
    std::unique_ptr<const Geometry> geom_;
    const oops::Variables vars_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_MODEL_MODEL_H_
