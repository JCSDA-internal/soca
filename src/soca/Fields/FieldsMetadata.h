/*
 * (C) Copyright 2024-2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

 */

#pragma once

#include <map>
#include <memory>
#include <string>

namespace soca {

// --------------------------------------------------------------------------------------

/// Holds all the user configurable meta data associated with a single field
/// NOTE: this is a subset of what exists on the Fortran side, and is being
/// expanded as it is needed during the Fortran -> C++ porting
struct FieldMetadata {
  std::string name;  // The internal (within soca) variable name
  std::string getvalName;  // The variable name used by UFO
  std::string getvalNameSurface;  // the variable name used by UFO for the surface
                                  // (if this is a 3D field)
};

// --------------------------------------------------------------------------------------

/// Metadata for the soca fields.
/// A yaml file is read in containing information about the various fields
class FieldsMetadata {
 public:
  /// Read metadata from a yaml file
  explicit FieldsMetadata(const std::string &);

  /// Get the metadata for a field
  /// @param name can either be the field's "name", "getvalName", or "getvalNameSurface"
  const FieldMetadata & operator[](const std::string & name) const;

 private:
  std::map<std::string, std::shared_ptr<FieldMetadata>> fieldMetadata_;
};

// --------------------------------------------------------------------------------------

}  // namespace soca
