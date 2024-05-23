/*
 * (C) Copyright 2024-2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

 */
#include <iostream>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Timer.h"

#include "soca/Fields/FieldsMetadata.h"

namespace soca {

FieldsMetadata::FieldsMetadata(const std::string & filename) {
  util::Timer timer("soca::FieldsMetadata", "FieldsMetadata");
  eckit::PathName filepathname = filename;
  eckit::YAMLConfiguration fullConfig(filepathname);

  // read in all the field metadata sections
  auto configs = fullConfig.getSubConfigurations("");
  for (auto config : configs) {
    auto f = std::make_shared<FieldMetadata>();

    f->name = config.getString("name");
    f->getvalName = config.getString("getval name", f->name);
    f->getvalNameSurface = config.getString("getval name surface", "");

    // check for duplicates
    if (fieldMetadata_.find(f->name) != fieldMetadata_.end() ||
        fieldMetadata_.find(f->getvalName) != fieldMetadata_.end() ||
        fieldMetadata_.find(f->getvalNameSurface) != fieldMetadata_.end() ) {
      util::abor1_cpp("Duplicate field metadata: "+f->name);
    }

    // insert
    fieldMetadata_[f->name] = f;
    fieldMetadata_[f->getvalName] = f;
    if (f->getvalNameSurface != "") {
      fieldMetadata_[f->getvalNameSurface] = f;
    }
  }
}

// --------------------------------------------------------------------------------------

const FieldMetadata & FieldsMetadata::operator[](const std::string &name) const {
  auto it = fieldMetadata_.find(name);
  if (it != fieldMetadata_.end()) {
    return *(it->second);
  } else {
    throw std::runtime_error("FieldMetadata not found for name: " + name);
  }
}

// --------------------------------------------------------------------------------------

}  // namespace soca
