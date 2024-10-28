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

/// Holds all the user configurable meta data associated with a single field.
/// This struct contains information about the internal variable name, the
/// variable name used by UFO, and the variable name used by UFO for the surface
/// (if this is a 3D field).
struct FieldMetadata {
  std::string name;  ///< The soca and JEDI name.
  std::string nameSurface;  ///< The variable name used by UFO for the surface
                            ///  (if this is a 3D field).
};

// --------------------------------------------------------------------------------------

/// Metadata for the soca fields.
/// A yaml file is read in containing information about the various fields.
class FieldsMetadata {
 public:
  /// Constructs a FieldsMetadata object and reads metadata from a yaml file.
  /// @param filepath The path to the yaml file.
  explicit FieldsMetadata(const std::string & filepath);

  /// Get the metadata for a field.
  /// @param name The field's "name", "getvalName", or "getvalNameSurface".
  /// @return A reference to the FieldMetadata object.
  const FieldMetadata & operator[](const std::string & name) const;

 private:
  /// Map of field name(s) to FieldMetadata objects.
  std::map<std::string, std::shared_ptr<FieldMetadata>> fieldMetadata_;
};

// --------------------------------------------------------------------------------------

}  // namespace soca
