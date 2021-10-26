/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_UFSM6C6_MODELUFSM6C6_H_
#define SOCA_MODEL_UFSM6C6_MODELUFSM6C6_H_

#include <memory>
#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "soca/Fortran.h"

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

  class ModelUFSm6c6:public oops::interface::ModelBase<Traits>,
              private util::ObjectCounter<ModelUFSm6c6>
  {
   public:
    static const std::string classname() {return "soca::ModelUFSm6c6";}

    ModelUFSm6c6(const Geometry &, const eckit::Configuration &);
    ~ModelUFSm6c6();

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
    std::unique_ptr<const Geometry> geom_;
    const oops::Variables vars_;
  };
  // -----------------------------------------------------------------------------

}  // namespace soca
#endif  // SOCA_MODEL_UFSM6C6_MODELUFSM6C6_H_
