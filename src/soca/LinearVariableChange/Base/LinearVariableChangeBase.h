/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/PolymorphicParameter.h"
#include "oops/util/Printable.h"

namespace soca {
  class Geometry;
  class Increment;
  class State;

// -----------------------------------------------------------------------------

class LinearVariableChangeParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(LinearVariableChangeParametersBase, oops::Parameters)
 public:
  oops::OptionalParameter<oops::Variables> inputVariables{"input variables", this};
  oops::OptionalParameter<oops::Variables> outputVariables{"output variables", this};
  oops::OptionalParameter<std::string> name{"linear variable change name", this};
};

// -----------------------------------------------------------------------------

class GenericLinearVariableChangeParameters :
                                     public LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericLinearVariableChangeParameters,
                           LinearVariableChangeParametersBase)
 public:
  oops::ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

class LinearVariableChangeBase : public util::Printable,
                                 private boost::noncopyable {
 public:
  LinearVariableChangeBase() {}
  virtual ~LinearVariableChangeBase() {}
  virtual void multiply(const Increment &, Increment &) const = 0;
  virtual void multiplyInverse(const Increment &, Increment &) const = 0;
  virtual void multiplyAD(const Increment &, Increment &) const = 0;
  virtual void multiplyInverseAD(const Increment &, Increment &) const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

class LinearVariableChangeFactory;

// -----------------------------------------------------------------------------

class LinearVariableChangeParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParametersWrapper, Parameters)
 public:
  oops::PolymorphicParameter<LinearVariableChangeParametersBase,
         LinearVariableChangeFactory>
  linearVariableChangeParameters{"linear variable change name",
                                 "default", this};
};

// -----------------------------------------------------------------------------

class LinearVariableChangeFactory {
 public:
  static LinearVariableChangeBase * create(const State & xbg, const State & xfg,
                                           const Geometry & geom,
                             const LinearVariableChangeParametersBase & params);

  static std::unique_ptr<LinearVariableChangeParametersBase>
         createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~LinearVariableChangeFactory() = default;

 protected:
  explicit LinearVariableChangeFactory(const std::string &name);

 private:
  virtual LinearVariableChangeBase * make(const State &, const State &,
                                          const Geometry &,
                                const LinearVariableChangeParametersBase &) = 0;

  virtual std::unique_ptr<LinearVariableChangeParametersBase>
          makeParameters() const = 0;

  static std::map < std::string, LinearVariableChangeFactory * > & getMakers() {
    static std::map < std::string, LinearVariableChangeFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class LinearVariableChangeMaker : public LinearVariableChangeFactory {
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T,
                     GenericLinearVariableChangeParameters>
    Parameters_;

  LinearVariableChangeBase * make(const State & xbg, const State & xfg,
                                  const Geometry & geom,
                  const LinearVariableChangeParametersBase & params) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(xbg, xfg, geom,
                oops::parametersOrConfiguration<oops::HasParameters_<T>::value>(
                   stronglyTypedParams));
  }

  std::unique_ptr<LinearVariableChangeParametersBase>
                  makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit LinearVariableChangeMaker(const std::string & name)
  : LinearVariableChangeFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace soca
