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
    // create a shared object that can be passed around
    auto f = std::make_shared<FieldMetadata>();

    // read in portion of the values (not all are read in on the C++ side, until
    // needed)
    f->name = config.getString("name");
    f->nameSurface = config.getString("name surface", "");

    // check for duplicates
    if (fieldMetadata_.count(f->name) > 0 ||
        fieldMetadata_.count(f->nameSurface) > 0) {
      util::abor1_cpp("Duplicate field metadata: " + f->name);
    }

    // insert into the maps, multiple copies are inserted for valid name
    fieldMetadata_[f->name] = f;
    if (f->nameSurface != "") {
      fieldMetadata_[f->nameSurface] = f;
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
