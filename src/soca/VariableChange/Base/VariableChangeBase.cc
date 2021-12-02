/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/VariableChange/Base/VariableChangeBase.h"

#include <vector>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace soca {

// -----------------------------------------------------------------------------

VariableChangeFactory::VariableChangeFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name
                        << " already registered in soca::VariableChangeFactory."
                       << std::endl;
    ABORT("Element already registered in soca::VariableChangeFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

VariableChangeBase * VariableChangeFactory::create(const Geometry & geom,
                                  const VariableChangeParametersBase & params) {
  oops::Log::trace() << "VariableChangeBase::create starting" << std::endl;

  const std::string &id = params.name.value().value();

  typename std::map<std::string, VariableChangeFactory*>::iterator
                                                    jloc = getMakers().find(id);

  if (jloc == getMakers().end()) {
    oops::Log::error() << id
              << " does not exist in soca::VariableChangeFactory." << std::endl;
    ABORT("Element does not exist in soca::VariableChangeFactory.");
  }

  VariableChangeBase * ptr = jloc->second->make(geom, params);
  oops::Log::trace() << "VariableChangeBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<VariableChangeParametersBase>
VariableChangeFactory::createParameters(const std::string &name) {
  typename std::map<std::string, VariableChangeFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name +
                              " does not exist in soca::VariableChangeFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace soca
