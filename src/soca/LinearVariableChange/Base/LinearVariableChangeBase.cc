/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/LinearVariableChange/Base/LinearVariableChangeBase.h"

#include <vector>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace soca {

// -------------------------------------------------------------------------------------------------

LinearVariableChangeFactory::LinearVariableChangeFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in soca::LinearVariableChangeFactory."
                       << std::endl;
    ABORT("Element already registered in soca::LinearVariableChangeFactory.");
  }
  getMakers()[name] = this;
}

// -------------------------------------------------------------------------------------------------

LinearVariableChangeBase * LinearVariableChangeFactory::create(const State & xbg,
     const State & xfg, const Geometry & geom, const LinearVariableChangeParametersBase & params) {
  oops::Log::trace() << "LinearVariableChangeBase::create starting" << std::endl;
  const std::string &id = params.name.value().value();
  typename std::map<std::string, LinearVariableChangeFactory*>::iterator
                                                                        jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in soca::LinearVariableChangeFactory."
                       << std::endl;
    ABORT("Element does not exist in soca::LinearVariableChangeFactory.");
  }
  LinearVariableChangeBase * ptr = jloc->second->make(xbg, xfg, geom, params);
  oops::Log::trace() << "LinearVariableChangeBase::create done" << std::endl;
  return ptr;
}

// -------------------------------------------------------------------------------------------------

std::unique_ptr<LinearVariableChangeParametersBase>
LinearVariableChangeFactory::createParameters(const std::string &name) {
  typename std::map<std::string, LinearVariableChangeFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in soca::LinearVariableChangeFactory");
  }
  return it->second->makeParameters();
}

// -------------------------------------------------------------------------------------------------

}  // namespace soca
